#include "ramcore/MateNames.h"
#include "rntuple/RAMNTupleRecord.h"

#include <ROOT/REntry.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleTypes.hxx>

#include <cstdint>
#include <string>
#include <utility>

void MateNameWriter::Move(RAMNTupleRecord &rec, ROOT::REntry &entry)
{
   if (!m_token)
      m_token = entry.GetToken(RAMNTupleRecord::kQnameRefField);
   std::uint32_t &ref = *entry.GetPtr<std::uint32_t>(*m_token);
   ref = 0;

   if (m_row - m_recentStart >= kMateNameWindow) {
      m_older = std::move(m_recent);
      m_recent.clear();
      m_recentStart = m_row;
   }

   // A name is referred to once: after the pair, a third record with it stores it again.
   auto *map = &m_recent;
   auto it = m_recent.find(rec.qname);
   if (it == m_recent.end()) {
      map = &m_older;
      it = m_older.find(rec.qname);
   }
   if (it != map->end()) {
      ref = static_cast<std::uint32_t>(m_row - it->second);
      map->erase(it);
      rec.qname.clear();
   } else {
      m_recent.emplace(rec.qname, m_row);
   }
   m_row++;
}

void MateNameWriter::Reset()
{
   m_recent.clear();
   m_older.clear();
   m_row = 0;
   m_recentStart = 0;
}

MateNameReader::MateNameReader(ROOT::RNTupleReader &reader)
{
   if (reader.GetDescriptor().FindFieldId(RAMNTupleRecord::kQnameRefField) == ROOT::kInvalidDescriptorId)
      return;
   m_ref.emplace(reader.GetView<std::uint32_t>(RAMNTupleRecord::kQnameRefField));
   m_names.emplace(reader.GetView<std::string>("record.qname"));
}

std::string MateNameReader::Get(const RAMNTupleRecord &rec, ROOT::NTupleSize_t row)
{
   const std::uint32_t ref = m_ref ? (*m_ref)(row) : 0;
   if (ref == 0)
      return rec.GetQNAME();
   return (*m_names)(row - ref);
}
