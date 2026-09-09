#include <benchmark/benchmark.h>
#include "benchmark_utils.h"
#include "ramcore/RAMNTupleView.h"
#include <string>
#include <sstream>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cstdint>
#include <memory>
#include <unistd.h>

// File and regions come from the environment so the numbers are reproducible
// on any machine. A hard-coded path is not a benchmark result.
//
//   RAMTOOLS_BENCH_RNTUPLE   RAM file written by samtoramntuple   (required)
//   RAMTOOLS_BENCH_REGIONS   comma-separated regions to query     (optional)

static std::vector<std::string> LoadRegions()
{
   const char *env = std::getenv("RAMTOOLS_BENCH_REGIONS");
   std::stringstream spec((env && *env) ? env
                                        : "chr1:1000000-1001000,chr1:1000000-2000000,chr1:1-50000000,chr21:1-48129895");
   std::vector<std::string> regions;
   std::string one{};
   while (std::getline(spec, one, ','))
      if (!one.empty())
         regions.push_back(one);
   return regions;
}

static const std::vector<std::string> kRegions = LoadRegions();

class RegionQueryFixture : public benchmark::Fixture {
public:
   void SetUp(const benchmark::State &state) override
   {
      region_idx_ = static_cast<int>(state.range(0));
      const char *value = std::getenv("RAMTOOLS_BENCH_RNTUPLE");
      rntuple_root_file_ = (value && *value) ? value : "";
   }

   void TearDown(const benchmark::State &) override {}

protected:
   int region_idx_;
   std::string rntuple_root_file_;

   // The query prints a summary per call. Reopening /dev/tty to undo the
   // redirect fails whenever the output is piped, so the descriptor is saved.
   void suppress_output()
   {
      fflush(stdout);
      m_saved_stdout = dup(STDOUT_FILENO);
      const std::unique_ptr<FILE, int (*)(FILE *)> null(fopen(NULL_DEVICE, "w"), fclose);
      dup2(fileno(null.get()), STDOUT_FILENO);
   }

   void restore_output() const
   {
      fflush(stdout);
      dup2(m_saved_stdout, STDOUT_FILENO);
      close(m_saved_stdout);
   }

   [[nodiscard]] const char *get_current_region() const { return kRegions[region_idx_].c_str(); }

private:
   int m_saved_stdout = -1;
};

// NOLINTNEXTLINE(misc-use-internal-linkage)
BENCHMARK_DEFINE_F(RegionQueryFixture, RNTuple)(benchmark::State &state)
{
   if (rntuple_root_file_.empty()) {
      state.SkipWithError("set RAMTOOLS_BENCH_RNTUPLE to a file written by samtoramntuple");
      return;
   }

   const char *region = get_current_region();
   int64_t total_reads_processed = 0;
   std::int64_t reads_in_this_run = 0;
   const RAMNTupleViewOpts opts{};

   for (auto _ : state) {
      suppress_output();
      reads_in_this_run = ramntupleview(rntuple_root_file_.c_str(), region, opts);
      restore_output();

      total_reads_processed += reads_in_this_run;
   }

   state.SetItemsProcessed(total_reads_processed);
   state.SetLabel(std::string(region) + ": " + std::to_string(reads_in_this_run) + " reads");
}

// One case per region. The cases used to be indices 0, 3, 6 and 9 into a
// fixed list of eighteen.
namespace {
BENCHMARK_REGISTER_F(RegionQueryFixture, RNTuple)
   ->DenseRange(0, static_cast<int>(kRegions.size()) - 1)
   ->Unit(benchmark::kSecond);
} // namespace

BENCHMARK_MAIN();
