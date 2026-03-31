#include <ROOT/RNTupleInspector.hxx>
#include <ROOT/RNTuple.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleDescriptor.hxx>
#include <TFile.h>
#include <iostream>
#include <iomanip>

void inspect_ram(const char *filename = "test/rntuple.root",
                 const char *ntupleName = "RAM")
{
   auto file = std::unique_ptr<TFile>(TFile::Open(filename));
   if (!file || !file->IsOpen()) {
      std::cerr << "Cannot open " << filename << "\n";
      return;
   }

   auto ntuple = file->Get<ROOT::RNTuple>(ntupleName);
   if (!ntuple) {
      std::cerr << "No RNTuple named '" << ntupleName << "' in file\n";
      std::cout << "\nKeys in file:\n";
      file->ls();
      return;
   }

   auto inspector = ROOT::Experimental::RNTupleInspector::Create(*ntuple);

   std::cout << "=== RNTuple: " << ntupleName << " ===\n";
   std::cout << "File:              " << filename << "\n";
   std::cout << "Entries:           " << inspector->GetDescriptor().GetNEntries() << "\n";
   std::cout << "Compressed size:   " << inspector->GetCompressedSize() << " bytes ("
             << std::fixed << std::setprecision(2)
             << inspector->GetCompressedSize() / 1024.0 / 1024.0 << " MB)\n";
   std::cout << "Uncompressed size: " << inspector->GetUncompressedSize() << " bytes ("
             << inspector->GetUncompressedSize() / 1024.0 / 1024.0 << " MB)\n";
   std::cout << "Compression ratio: " << std::setprecision(2)
             << inspector->GetCompressionFactor() << "x\n\n";

   // Column type summary — this is the key output
   std::cout << "=== Column Type Summary ===\n";
   inspector->PrintColumnTypeInfo();
   std::cout << "\n";

   // Field + column structure from descriptor
   std::cout << "=== Field Structure ===\n";
   const auto &descriptor = inspector->GetDescriptor();
   for (const auto &fieldDesc : descriptor.GetTopLevelFields()) {
      std::cout << fieldDesc.GetFieldName()
                << " [" << fieldDesc.GetTypeName() << "]\n";
      for (const auto &subField : descriptor.GetFieldIterable(fieldDesc.GetId())) {
         std::cout << "  └─ " << subField.GetFieldName()
                   << " [" << subField.GetTypeName() << "]\n";
      }
   }

   // All columns
   std::cout << "\n=== All Columns ===\n";
   std::cout << std::left
             << std::setw(8)  << "ColID"
             << std::setw(10) << "FieldID"
             << std::setw(30) << "FieldName"
             << std::setw(20) << "ColType"
             << "\n";
   std::cout << std::string(68, '-') << "\n";

   for (const auto &fieldDesc : descriptor.GetFieldIterable(descriptor.GetFieldZeroId())) {
      for (const auto &colDesc : descriptor.GetColumnIterable(fieldDesc.GetId())) {
         std::cout << std::left
                   << std::setw(8)  << colDesc.GetPhysicalId()
                   << std::setw(10) << colDesc.GetFieldId()
                   << std::setw(30) << fieldDesc.GetFieldName()
                   << std::setw(20) << static_cast<int>(colDesc.GetType())
                   << "\n";
      }
   }

   // METADATA and INDEX
   const char *extras[] = {"METADATA", "INDEX"};
   for (const char *name : extras) {
      auto extra = file->Get<ROOT::RNTuple>(name);
      if (!extra) continue;

      auto extraInsp = ROOT::Experimental::RNTupleInspector::Create(*extra);
      std::cout << "\n=== RNTuple: " << name << " ===\n";
      std::cout << "Entries:           " << extraInsp->GetDescriptor().GetNEntries() << "\n";
      std::cout << "Compressed size:   " << extraInsp->GetCompressedSize() << " bytes\n";
      std::cout << "Uncompressed size: " << extraInsp->GetUncompressedSize() << " bytes\n";
   }

   std::cout << "\n=== File size on disk: "
             << file->GetSize() / 1024.0 / 1024.0 << " MB ===\n";
}