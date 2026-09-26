#ifndef RAMCORE_QUALITYBLOCKS_H
#define RAMCORE_QUALITYBLOCKS_H

// Quality scores compressed with fqzcomp in blocks of records. The last record
// of a block carries it in "qualblock", and the block's SEQ (SeqBlocks.h) in
// "seqblock"; METADATA lists the row each block ends on.

#include "ramcore/BlockReadAhead.h"
#include "ramcore/SeqBlocks.h"
#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/REntry.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleTypes.hxx>
#include <ROOT/RNTupleView.hxx>

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <vector>

constexpr std::size_t kQualityBlockRecords = 25000;

/// Holds one record back so the record that ends a block can carry the block's
/// qualities and SEQ. Only Phred+33 qualities in SAM's '!'..'~' go into a
/// block, and only SEQ the model codes; the rest stay in the record.
class QualityBlockWriter {
public:
   using FillFn = std::function<void(ROOT::REntry &)>;

   QualityBlockWriter(std::unique_ptr<ROOT::REntry> first, std::unique_ptr<ROOT::REntry> second, FillFn fill,
                      std::size_t block_records = kQualityBlockRecords);

   RAMNTupleRecord &Record() { return *m_records.at(m_current); }
   /// The entry that holds Record(), for other fields of the same row.
   ROOT::REntry &Entry() { return *m_entries.at(m_current); }
   void Add();
   /// Ends the current block; needed before the written rows are flushed.
   void Finish();
   /// Last row of each block finished since the previous call, counted from the
   /// first row after it.
   std::vector<uint64_t> TakeBlockEnds();

private:
   std::array<std::unique_ptr<ROOT::REntry>, 2> m_entries;
   std::array<RAMNTupleRecord *, 2> m_records{};
   std::array<std::vector<std::uint8_t> *, 2> m_blobs{};
   std::array<std::vector<std::uint8_t> *, 2> m_seqBlobs{};
   FillFn m_fill;
   std::size_t m_blockRecords;
   std::size_t m_current = 0;
   bool m_pending = false;
   std::size_t m_rows = 0;
   uint64_t m_rowsSinceTake = 0;
   std::string m_quals;
   std::vector<uint32_t> m_lengths;
   std::vector<uint32_t> m_flags;
   std::vector<uint64_t> m_blockEnds;
   SeqBlockEncoder m_seqs;
};

/// Returns QUAL as SAM text, decoding one block at a time and decoding ahead on
/// other threads during a scan. Needs the metadata loaded by OpenRAMFile().
class QualityBlockReader {
public:
   explicit QualityBlockReader(ROOT::RNTupleReader &reader);

   std::string Get(const RAMNTupleRecord &rec, ROOT::NTupleSize_t row);

private:
   struct Packed {
      std::vector<std::size_t> slots; ///< Per row: index of its quality in the block.
      std::size_t records = 0;
      std::vector<char> bytes;
   };
   struct Block {
      std::vector<std::size_t> slots;
      std::string quals;
      std::vector<int> lengths;
      std::vector<std::size_t> offsets;
   };

   Packed Read(std::size_t block);
   static Block Decode(Packed packed);

   ROOT::RNTupleView<uint32_t> m_flagsView;
   std::optional<ROOT::RNTupleView<std::vector<std::uint8_t>>> m_view;
   BlockReadAhead<Packed, Block> m_blocks{&Decode};
};

#endif // RAMCORE_QUALITYBLOCKS_H
