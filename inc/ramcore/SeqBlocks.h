#ifndef RAMCORE_SEQBLOCKS_H
#define RAMCORE_SEQBLOCKS_H

// SEQ compressed per record block with a context model: each base is predicted
// from the bases earlier reads of the block put at the same reference
// position, else from an earlier copy of the preceding bases, else from the
// preceding 11 bases. The last record of a block carries it in "seqblock",
// next to its "qualblock".

#include "ramcore/BlockReadAhead.h"
#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleTypes.hxx>
#include <ROOT/RNTupleView.hxx>

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

/// What the model needs of a read besides its bases.
struct SeqBlockRead {
   int32_t refid = -1;
   int32_t pos = -1;
   std::vector<uint32_t> cigar;
};

/// Collects the SEQ of a block's records and codes them. Takes reads whose SEQ
/// is only A, C, G, T and N, all of which come back exactly; other reads keep
/// their SEQ in the record.
class SeqBlockEncoder {
public:
   /// Moves rec's SEQ into the block and sets kSeqInBlock if it can be coded,
   /// else clears kSeqInBlock and leaves the SEQ in the record.
   void Add(RAMNTupleRecord &rec);
   /// The coded block, empty if no read went in; starts a new block.
   std::vector<std::uint8_t> Finish();

private:
   std::vector<SeqBlockRead> m_reads;
   std::vector<std::string> m_seqs;
};

/// Returns SEQ as SAM text. Needs the metadata loaded by OpenRAMFile().
class SeqBlockReader {
public:
   explicit SeqBlockReader(ROOT::RNTupleReader &reader);

   std::string Get(const RAMNTupleRecord &rec, ROOT::NTupleSize_t row);

private:
   struct Packed {
      std::vector<std::size_t> slots; ///< Per row: index of its SEQ in the block.
      std::vector<SeqBlockRead> reads;
      std::vector<std::uint8_t> bytes;
   };
   struct Block {
      std::vector<std::size_t> slots;
      std::vector<std::string> seqs;
   };

   Packed Read(std::size_t block);
   static Block Decode(Packed packed);

   ROOT::RNTupleView<uint32_t> m_flagsView;
   ROOT::RNTupleView<int32_t> m_refidView;
   ROOT::RNTupleView<int32_t> m_posView;
   ROOT::RNTupleView<std::vector<uint32_t>> m_cigarView;
   std::optional<ROOT::RNTupleView<std::vector<std::uint8_t>>> m_view;
   BlockReadAhead<Packed, Block> m_blocks{&Decode};
};

namespace SeqBlockCodec {
/// True if the model can code \p seq: not empty, not "*", only A, C, G, T and N.
bool Codable(const std::string &seq);
std::vector<std::uint8_t> Encode(const std::vector<SeqBlockRead> &reads, const std::vector<std::string> &seqs);
/// Throws if \p bytes is not a block coded for \p reads.
std::vector<std::string> Decode(const std::vector<std::uint8_t> &bytes, const std::vector<SeqBlockRead> &reads);
} // namespace SeqBlockCodec

#endif // RAMCORE_SEQBLOCKS_H
