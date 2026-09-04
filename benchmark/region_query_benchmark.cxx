#include <benchmark/benchmark.h>
#include "benchmark_utils.h"
#include "ramcore/RAMNTupleView.h"
#include <string>
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <cstdint>

// File comes from the environment so the numbers are reproducible on any
// machine. A hard-coded path is not a benchmark result.
//
//   RAMTOOLS_BENCH_RNTUPLE   RAM file written by samtoramntuple   (required)
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

   static const std::vector<std::string> regions_;

   void suppress_output() { freopen(NULL_DEVICE, "w", stdout); }
   void restore_output() { freopen("/dev/tty", "w", stdout); }

   const char *get_current_region() const { return regions_[region_idx_ % regions_.size()].c_str(); }
};

const std::vector<std::string> RegionQueryFixture::regions_ = {"chr1:1000000-1001000",
                                                               "chr2:5000000-5010000",
                                                               "chrX:100000-150000",
                                                               "chr1:1000000-2000000",
                                                               "chr5:10000000-15000000",
                                                               "chr10:50000000-60000000",
                                                               "chr1:1-50000000",
                                                               "chr2:1-100000000",
                                                               "chr7:50000000-150000000",
                                                               "chr21:1-48129895",
                                                               "chrM:1-16571",
                                                               "chrY:2600000-2700000",
                                                               "GL000227.1:1-100000",
                                                               "chr1:1-1000",
                                                               "chr1:249250621-249250621",
                                                               "chr22:51304566-51304566",
                                                               "chr17:41196312-41277500",
                                                               "chr13:32889611-32973805"};

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
   state.counters["region_idx"] = region_idx_;
   state.SetLabel(std::to_string(reads_in_this_run) + " reads");
}

BENCHMARK_REGISTER_F(RegionQueryFixture, RNTuple)->Args({0})->Args({3})->Args({6})->Args({9})->Unit(benchmark::kSecond);

BENCHMARK_MAIN();
