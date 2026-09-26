#include "ramcore/SeqBlocks.h"
#include "ramcore/CigarOps.h"
#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleTypes.hxx>

#include <algorithm>
#include <array>
#include <cstddef>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr std::uint8_t kFormat = 1;
constexpr uint16_t kHalf = 1 << 15;
constexpr int kShiftFast = 4;
constexpr int kShiftSlow = 5;
constexpr int kOrder = 11;    // preceding bases the fallback context holds
constexpr int kMatchMin = 12; // bases that must repeat before a match is followed
// Tables sized for a block of 25,000 short reads; larger ones gain under 1% and decode slower.
constexpr int kHashBits = 18;
constexpr int kPileBits = 18;
// Guards the decoder against a corrupt length; far above any real read.
constexpr uint32_t kMaxLength = 1U << 28;

/// Binary range coder (LZMA's): Code(p, shift, bit) codes bit with a
/// probability of p / 65536 for a 0 and moves p towards the bit coded.
class RangeEncoder {
public:
   static constexpr bool kEncoding = true;

   int Code(uint16_t &p, int shift, int bit)
   {
      const uint32_t bound = (m_range >> 16) * p;
      if (bit == 0) {
         m_range = bound;
         p += (65536 - p) >> shift;
      } else {
         m_low += bound;
         m_range -= bound;
         p -= p >> shift;
      }
      while (m_range < (1U << 24)) {
         m_range <<= 8;
         ShiftLow();
      }
      return bit;
   }

   void Finish(std::vector<std::uint8_t> &out)
   {
      for (int i = 0; i < 5; i++)
         ShiftLow();
      out.insert(out.end(), m_out.begin(), m_out.end());
   }

private:
   void ShiftLow()
   {
      if (static_cast<uint32_t>(m_low) < 0xFF000000U || (m_low >> 32) != 0) {
         const auto carry = static_cast<std::uint8_t>(m_low >> 32);
         std::uint8_t pending = m_cache;
         for (; m_cacheSize > 0; m_cacheSize--) {
            m_out.push_back(static_cast<std::uint8_t>(pending + carry));
            pending = 0xFF;
         }
         m_cache = static_cast<std::uint8_t>(m_low >> 24);
      }
      m_cacheSize++;
      m_low = (m_low & 0x00FFFFFF) << 8;
   }

   uint64_t m_low = 0;
   uint32_t m_range = 0xFFFFFFFF;
   std::uint8_t m_cache = 0;
   uint64_t m_cacheSize = 1;
   std::vector<std::uint8_t> m_out;
};

class RangeDecoder {
public:
   static constexpr bool kEncoding = false;

   /// Decodes \p bytes from \p start on.
   RangeDecoder(const std::vector<std::uint8_t> &bytes, std::size_t start) : m_bytes(&bytes), m_next(start)
   {
      for (int i = 0; i < 5; i++)
         m_code = (m_code << 8) | Next();
   }

   int Code(uint16_t &p, int shift, int /*bit*/)
   {
      const uint32_t bound = (m_range >> 16) * p;
      int bit = 0;
      if (m_code < bound) {
         m_range = bound;
         p += (65536 - p) >> shift;
      } else {
         m_code -= bound;
         m_range -= bound;
         p -= p >> shift;
         bit = 1;
      }
      while (m_range < (1U << 24)) {
         m_range <<= 8;
         m_code = (m_code << 8) | Next();
      }
      return bit;
   }

private:
   std::uint8_t Next() { return m_next < m_bytes->size() ? (*m_bytes)[m_next++] : 0; }

   const std::vector<std::uint8_t> *m_bytes;
   std::size_t m_next;
   uint32_t m_range = 0xFFFFFFFF;
   uint32_t m_code = 0;
};

constexpr std::array<char, 4> kBases{'A', 'C', 'G', 'T'};

int BaseCode(char c)
{
   switch (c) {
   case 'A': return 0;
   case 'C': return 1;
   case 'G': return 2;
   case 'T': return 3;
   default: return -1;
   }
}

/// Bases of the read the CIGAR accounts for; 0 without a CIGAR.
uint32_t QueryLength(const std::vector<uint32_t> &cigar)
{
   uint32_t length = 0;
   for (const uint32_t c : cigar) {
      const uint32_t op = c & 0xF;
      if (op == RAM_CIGAR_M || op == RAM_CIGAR_I || op == RAM_CIGAR_S || op == RAM_CIGAR_EQUAL || op == RAM_CIGAR_X)
         length += c >> 4;
   }
   return length;
}

/// Reference position of each base, -1 where it has none (insertions, clips, unplaced reads).
void RefPositions(const SeqBlockRead &read, std::size_t length, std::vector<int64_t> &out)
{
   out.assign(length, -1);
   if (read.refid < 0 || read.pos < 0)
      return;
   int64_t ref = read.pos;
   std::size_t q = 0;
   for (const uint32_t c : read.cigar) {
      const uint32_t op = c & 0xF;
      const uint32_t n = c >> 4;
      if (op == RAM_CIGAR_M || op == RAM_CIGAR_EQUAL || op == RAM_CIGAR_X) {
         for (uint32_t i = 0; i < n && q < length; i++)
            out[q++] = ref++;
      } else if (op == RAM_CIGAR_I || op == RAM_CIGAR_S) {
         q += n;
      } else if (op == RAM_CIGAR_D || op == RAM_CIGAR_N) {
         ref += n;
      }
   }
}

/// Model state for one block; every block starts from a fresh one.
struct Model {
   // Three probabilities per context: the base's first bit, then its second after a 0 or a 1.
   /// Context: the base the pileup predicts, how often it was seen (up to 3), and whether the previous base missed.
   std::vector<uint16_t> pile = std::vector<uint16_t>(3 * 4 * 4 * 2, kHalf);
   /// Context: the base the match predicts and the match length in steps of 4 bases (up to 3).
   std::vector<uint16_t> match = std::vector<uint16_t>(3 * 4 * 4, kHalf);
   /// Context: a hash of the preceding kOrder bases.
   std::vector<uint16_t> order = std::vector<uint16_t>(std::size_t(3) << kHashBits, kHalf);

   /// Last base coded at a reference position, and how many reads in a row agreed on it.
   struct Cell {
      uint64_t key = ~0ULL;
      std::uint8_t base = 0;
      std::uint8_t seen = 0;
   };
   std::vector<Cell> pileup = std::vector<Cell>(std::size_t(1) << kPileBits);
   /// Per hash of kMatchMin bases: the position in `bases` that last followed them.
   std::vector<uint32_t> follows = std::vector<uint32_t>(std::size_t(1) << kHashBits, 0);
   /// Every base coded so far; position 0 means "none" in `follows`.
   std::vector<std::uint8_t> bases{0};
   std::vector<int64_t> refs;

   uint16_t sameLength = kHalf;
   uint32_t lastLength = 0;
   uint16_t hasN = kHalf;
   std::array<uint16_t, 2> isN{kHalf, kHalf}; ///< By whether the previous base was N.

   static uint64_t Key(int32_t refid, int64_t pos) { return (uint64_t(uint32_t(refid)) << 40) ^ uint64_t(pos); }
   static std::size_t Hash(uint64_t h, int bits) { return std::size_t((h * 0x9E3779B97F4A7C15ULL) >> (64 - bits)); }
   Cell &PileupCell(uint64_t key) { return pileup[key & ((uint64_t(1) << kPileBits) - 1)]; }
};

/// Codes the read's length as "what its CIGAR says" (the previous read's length
/// without a CIGAR) or 32 plain bits. Returns the length.
template <typename Coder>
uint32_t CodeLength(Coder &rc, Model &m, const SeqBlockRead &read, uint32_t length)
{
   const uint32_t query = QueryLength(read.cigar);
   const uint32_t predicted = query > 0 ? query : m.lastLength;
   if (rc.Code(m.sameLength, kShiftSlow, length == predicted)) {
      length = predicted;
   } else {
      uint32_t coded = 0;
      for (int k = 31; k >= 0; k--) {
         uint16_t half = kHalf;
         coded |= static_cast<uint32_t>(rc.Code(half, kShiftSlow, static_cast<int>((length >> k) & 1))) << k;
      }
      length = coded;
   }
   if (length == 0 || length > kMaxLength)
      throw std::runtime_error("SEQ block holds a read of length " + std::to_string(length));
   m.lastLength = length;
   return length;
}

/// Codes seq, which has the read's length; decoding overwrites it.
template <typename Coder>
void CodeBases(Coder &rc, Model &m, const SeqBlockRead &read, std::string &seq)
{
   const bool hasN = rc.Code(m.hasN, kShiftSlow, seq.find('N') != std::string::npos) != 0;
   RefPositions(read, seq.size(), m.refs);
   uint64_t hist = 0;
   std::size_t matchPtr = 0;
   std::size_t candidate = 0; // position that followed the last kMatchMin bases before
   int matchLen = 0;
   int lastMiss = 0;
   int prevN = 0;
   for (std::size_t i = 0; i < seq.size(); i++) {
      if (hasN) {
         prevN = rc.Code(m.isN.at(static_cast<std::size_t>(prevN)), kShiftFast, seq[i] == 'N');
         if (prevN) {
            seq[i] = 'N';
            matchPtr = 0;
            candidate = 0;
            continue;
         }
      }
      if (!matchPtr && candidate) {
         matchPtr = candidate;
         matchLen = 0;
      }
      candidate = 0;

      const int64_t ref = m.refs[i];
      const uint64_t key = ref >= 0 ? Model::Key(read.refid, ref) : 0;
      Model::Cell *cell = nullptr;
      if (ref >= 0 && m.PileupCell(key).key == key)
         cell = &m.PileupCell(key);

      std::vector<uint16_t> *probs = &m.order;
      std::size_t at = 0;
      int shift = kShiftSlow;
      int pred = -1;
      if (cell) {
         pred = cell->base;
         probs = &m.pile;
         at = 3 * static_cast<std::size_t>((pred * 4 + std::min<int>(cell->seen, 3)) * 2 + lastMiss);
      } else if (matchPtr) {
         pred = m.bases[matchPtr];
         probs = &m.match;
         at = 3 * static_cast<std::size_t>(pred * 4 + std::min(matchLen / 4, 3));
      } else {
         at = 3 * Model::Hash(hist & ((1ULL << (2 * kOrder)) - 1), kHashBits);
         shift = kShiftFast;
      }

      const int base = Coder::kEncoding ? BaseCode(seq[i]) : 0;
      const int hi = rc.Code((*probs)[at], shift, base >> 1);
      const int b = (hi << 1) | rc.Code((*probs)[at + 1 + static_cast<std::size_t>(hi)], shift, base & 1);
      seq[i] = kBases.at(static_cast<std::size_t>(b));

      lastMiss = pred >= 0 && b != pred;
      if (matchPtr) {
         if (m.bases[matchPtr] == b) {
            matchLen++;
            matchPtr++;
         } else {
            matchPtr = 0;
         }
      }
      if (ref >= 0) {
         Model::Cell &c = m.PileupCell(key);
         if (c.key == key && c.base == b) {
            if (c.seen < 255)
               c.seen++;
         } else {
            c = {key, static_cast<std::uint8_t>(b), 0};
         }
      }
      hist = (hist << 2) | static_cast<uint64_t>(b);
      m.bases.push_back(static_cast<std::uint8_t>(b));
      if (i + 1 >= kMatchMin) {
         // Look the bases up before recording where they are now, or the lookup finds itself.
         uint32_t &slot = m.follows[Model::Hash(hist & ((1ULL << (2 * kMatchMin)) - 1), kHashBits)];
         candidate = slot;
         slot = static_cast<uint32_t>(m.bases.size());
      }
   }
}

} // namespace

namespace SeqBlockCodec {

bool Codable(const std::string &seq)
{
   return !seq.empty() && seq != "*" && seq.size() <= kMaxLength &&
          std::all_of(seq.begin(), seq.end(), [](char c) { return c == 'N' || BaseCode(c) >= 0; });
}

std::vector<std::uint8_t> Encode(const std::vector<SeqBlockRead> &reads, const std::vector<std::string> &seqs)
{
   if (reads.size() != seqs.size())
      throw std::invalid_argument("SEQ block: reads and sequences differ in number");
   RangeEncoder rc;
   Model m;
   std::string seq;
   for (std::size_t i = 0; i < reads.size(); i++) {
      if (!Codable(seqs[i]))
         throw std::invalid_argument("SEQ block: cannot code \"" + seqs[i] + "\"");
      seq = seqs[i];
      CodeLength(rc, m, reads[i], static_cast<uint32_t>(seq.size()));
      CodeBases(rc, m, reads[i], seq);
   }
   std::vector<std::uint8_t> out{kFormat};
   rc.Finish(out);
   return out;
}

std::vector<std::string> Decode(const std::vector<std::uint8_t> &bytes, const std::vector<SeqBlockRead> &reads)
{
   if (bytes.empty() || bytes[0] != kFormat)
      throw std::runtime_error("SEQ block has an unknown format");
   RangeDecoder rc(bytes, /*start=*/1);
   Model m;
   std::vector<std::string> seqs(reads.size());
   for (std::size_t i = 0; i < reads.size(); i++) {
      seqs[i].assign(CodeLength(rc, m, reads[i], /*length=*/0), 'A');
      CodeBases(rc, m, reads[i], seqs[i]);
   }
   return seqs;
}

} // namespace SeqBlockCodec

void SeqBlockEncoder::Add(RAMNTupleRecord &rec)
{
   if (rec.TestBit(RAMNTupleRecord::kSeqRaw) && SeqBlockCodec::Codable(rec.seq)) {
      m_reads.push_back({rec.refid, rec.pos, rec.cigar});
      m_seqs.push_back(std::move(rec.seq));
      rec.seq.clear();
      rec.SetBit(RAMNTupleRecord::kSeqInBlock);
   } else {
      rec.compression_flags &= ~static_cast<uint32_t>(RAMNTupleRecord::kSeqInBlock);
   }
}

std::vector<std::uint8_t> SeqBlockEncoder::Finish()
{
   std::vector<std::uint8_t> bytes;
   if (!m_seqs.empty())
      bytes = SeqBlockCodec::Encode(m_reads, m_seqs);
   m_reads.clear();
   m_seqs.clear();
   return bytes;
}

SeqBlockReader::SeqBlockReader(ROOT::RNTupleReader &reader)
   : m_flagsView(reader.GetView<uint32_t>("record.compression_flags")),
     m_refidView(reader.GetView<int32_t>("record.refid")),
     m_posView(reader.GetView<int32_t>("record.pos")),
     m_cigarView(reader.GetView<std::vector<uint32_t>>("record.cigar"))
{
   if (reader.GetDescriptor().FindFieldId(RAMNTupleRecord::kSeqBlockField) != ROOT::kInvalidDescriptorId)
      m_view.emplace(reader.GetView<std::vector<std::uint8_t>>(RAMNTupleRecord::kSeqBlockField));
}

std::string SeqBlockReader::Get(const RAMNTupleRecord &rec, ROOT::NTupleSize_t row)
{
   if (!rec.TestBit(RAMNTupleRecord::kSeqInBlock))
      return rec.GetSEQ();

   uint64_t first = 0;
   const Block *block = m_view ? m_blocks.Find(row, first, [this](std::size_t b) { return Read(b); }) : nullptr;
   if (!block)
      throw std::runtime_error("record " + std::to_string(row) + " has no SEQ block");
   return block->seqs[block->slots[static_cast<std::size_t>(row - first)]];
}

SeqBlockReader::Packed SeqBlockReader::Read(std::size_t block)
{
   const auto &ends = RAMNTupleRecord::GetQualBlockEnds();
   const uint64_t first = block == 0 ? 0 : ends[block - 1] + 1;
   Packed p;
   p.slots.assign(static_cast<std::size_t>(ends[block] - first + 1), 0);
   for (std::size_t i = 0; i < p.slots.size(); i++) {
      const uint64_t row = first + i;
      if (m_flagsView(row) & RAMNTupleRecord::kSeqInBlock) {
         p.slots[i] = p.reads.size();
         p.reads.push_back({m_refidView(row), m_posView(row), m_cigarView(row)});
      }
   }
   if (!p.reads.empty())
      p.bytes = (*m_view)(ends[block]);
   return p;
}

SeqBlockReader::Block SeqBlockReader::Decode(Packed packed)
{
   Block b;
   b.slots = std::move(packed.slots);
   if (!packed.reads.empty())
      b.seqs = SeqBlockCodec::Decode(packed.bytes, packed.reads);
   return b;
}
