#ifndef RAMCORE_BLOCKREADAHEAD_H
#define RAMCORE_BLOCKREADAHEAD_H

// Decoded record blocks for a reader: the block the rows being read fall in,
// and during a scan the next ones, decoded on other threads.

#include "rntuple/RAMNTupleRecord.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <future>
#include <map>
#include <thread>
#include <utility>

/// Packed is what a block's rows hold, read on the caller's thread because views
/// are not thread-safe; Decode turns it into a Block and may run on any thread.
/// Needs the block ends loaded by OpenRAMFile().
template <typename Packed, typename Block>
class BlockReadAhead {
public:
   using DecodeFn = Block (*)(Packed);

   explicit BlockReadAhead(DecodeFn decode)
      : m_decode(decode), m_readAhead(std::max(1U, std::thread::hardware_concurrency()))
   {
   }

   /// The block that holds \p row, or nullptr if none does; \p first_row is set
   /// to the block's first row. \p read(block) returns the block's Packed.
   template <typename ReadFn>
   const Block *Find(uint64_t row, uint64_t &first_row, ReadFn &&read)
   {
      const auto &ends = RAMNTupleRecord::GetQualBlockEnds();
      const auto it = std::lower_bound(ends.begin(), ends.end(), row);
      if (it == ends.end())
         return nullptr;
      const auto block = static_cast<std::size_t>(it - ends.begin());
      if (block != m_block)
         Load(block, read);
      first_row = m_firstRow;
      return &m_current;
   }

private:
   static constexpr std::size_t kNoBlock = static_cast<std::size_t>(-1);

   template <typename ReadFn>
   void Load(std::size_t block, ReadFn &read)
   {
      const auto &ends = RAMNTupleRecord::GetQualBlockEnds();
      // Read ahead one block more for each consecutive block, so random lookups start nothing.
      m_run = (m_block != kNoBlock && block == m_block + 1) ? m_run + 1 : 0;

      auto ahead = m_ahead.find(block);
      if (ahead != m_ahead.end()) {
         m_current = ahead->second.get();
         m_ahead.erase(ahead);
      } else {
         m_current = m_decode(read(block));
      }
      m_block = block;
      m_firstRow = block == 0 ? 0 : ends[block - 1] + 1;

      const std::size_t depth = std::min(m_run, m_readAhead);
      for (std::size_t next = block + 1; next < ends.size() && next <= block + depth; next++) {
         if (m_ahead.count(next) == 0)
            m_ahead.emplace(next, std::async(std::launch::async, m_decode, read(next)));
      }
   }

   DecodeFn m_decode;
   std::size_t m_readAhead;
   std::size_t m_block = kNoBlock;
   std::size_t m_run = 0; ///< Consecutive blocks read so far.
   uint64_t m_firstRow = 0;
   Block m_current;
   std::map<std::size_t, std::future<Block>> m_ahead;
};

#endif // RAMCORE_BLOCKREADAHEAD_H
