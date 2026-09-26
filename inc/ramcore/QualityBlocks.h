#ifndef RAMCORE_QUALITYBLOCKS_H
#define RAMCORE_QUALITYBLOCKS_H

// Quality scores compressed with fqzcomp, the CRAM 3.1 quality codec, in blocks
// of records. The last record of each block carries the compressed block in the
// "qualblock" field; METADATA lists the row each block ends on.

#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/REntry.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleView.hxx>

#include <cstddef>
#include <cstdint>
#include <functional>
#include <future>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

/// Records per quality block. Larger blocks compress better; a region dump
/// decodes whole blocks.
constexpr std::size_t kQualityBlockRecords = 25000;

/// Fills records through two entries, holding the last one back, so that the
/// record that ends a block can carry the block. Qualities stored as Phred+33
/// text go into the block; "*", any other policy, or text outside SAM's
/// '!'..'~' stays in the record.
class QualityBlockWriter {
public:
   using FillFn = std::function<void(ROOT::REntry &)>;

   QualityBlockWriter(std::unique_ptr<ROOT::REntry> first, std::unique_ptr<ROOT::REntry> second, FillFn fill,
                      std::size_t block_records = kQualityBlockRecords);

   /// The record to set for the next input record.
   RAMNTupleRecord &Record() { return *fRecords[fCurrent]; }
   /// Takes the record set through Record().
   void Add();
   /// Ends the current block. Call it at the end of the input and wherever the
   /// rows written so far must be complete, such as before a cluster is flushed.
   void Finish();
   /// Records in each block finished since the last TakeBlockSizes(), in order.
   std::vector<uint32_t> TakeBlockSizes();

private:
   std::unique_ptr<ROOT::REntry> fEntries[2];
   RAMNTupleRecord *fRecords[2];
   std::vector<std::uint8_t> *fBlobs[2];
   FillFn fFill;
   std::size_t fBlockRecords;
   int fCurrent = 0;
   bool fPending = false;
   std::size_t fRows = 0;
   std::string fQuals;
   std::vector<uint32_t> fLengths;
   std::vector<uint32_t> fFlags;
   std::vector<uint32_t> fBlockSizes;
};

/// Gives back QUAL as SAM text for records read from a RAM file, decoding one
/// block at a time. Needs RAMNTupleRecord::OpenRAMFile() to have loaded the
/// file's metadata. Reading consecutive blocks decodes the next ones on other
/// threads.
class QualityBlockReader {
public:
   explicit QualityBlockReader(ROOT::RNTupleReader &reader);

   std::string Get(const RAMNTupleRecord &rec, ROOT::NTupleSize_t row);

private:
   struct Packed {
      std::vector<std::size_t> slots; ///< Per row of the block: its string, if in the block.
      std::size_t records = 0;        ///< Qualities in the block.
      std::vector<char> bytes;
   };
   struct Block {
      std::vector<std::size_t> slots;
      std::string quals;
      std::vector<int> lengths;
      std::vector<std::size_t> offsets;
   };
   static constexpr std::size_t kNoBlock = static_cast<std::size_t>(-1);

   /// Reads a block's flags and bytes; views are not thread-safe, so the caller's thread does this.
   Packed Read(std::size_t block);
   /// Decodes a block; safe on any thread.
   static Block Decode(Packed packed);
   void Load(std::size_t block);

   ROOT::RNTupleView<uint32_t> fFlagsView;
   std::optional<ROOT::RNTupleView<std::vector<std::uint8_t>>> fView;
   std::size_t fReadAhead;
   std::size_t fBlock = kNoBlock;
   std::size_t fRun = 0; ///< Consecutive blocks read so far.
   uint64_t fFirstRow = 0;
   Block fCurrent;
   std::map<std::size_t, std::future<Block>> fAhead;
};

#endif // RAMCORE_QUALITYBLOCKS_H
