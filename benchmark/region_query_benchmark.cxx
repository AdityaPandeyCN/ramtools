#include <benchmark/benchmark.h>
#include "benchmark_utils.h"
#include "ramcore/RAMNTupleView.h"
#include "ramcore/SamToTTree.h"
#include "ramcore/SamToNTuple.h"
#include <string>
#include <cstdio>
#include <cstdlib>
#include <unistd.h>
#include <vector>
#include <Rtypes.h>

Long64_t ramview(const char *file, const char *query, bool cache = true, bool perfstats = false,
                 const char *perfstatsfilename = "perf.root");

// The files and the regions come from the environment, because the answer this
// benchmark gives depends entirely on which alignment file it was pointed at --
// a hard-coded path to one machine's disk makes the numbers unreproducible, and
// unreproducible numbers are not a benchmark result.
//
//   RAMTOOLS_BENCH_RNTUPLE   RAM file written by samtoramntuple   (required)
//   RAMTOOLS_BENCH_TTREE     RAM file written by samtoram         (optional)
//   RAMTOOLS_BENCH_REGIONS   comma-separated regions to query     (optional)
class RegionQueryFixture : public benchmark::Fixture {
public:
   void SetUp(const benchmark::State &state) override
   {
      region_idx_ = static_cast<int>(state.range(0));

      ttree_root_file_ = EnvOr("RAMTOOLS_BENCH_TTREE", "");
      rntuple_root_file_ = EnvOr("RAMTOOLS_BENCH_RNTUPLE", "");
   }

   void TearDown(const benchmark::State &) override {}

protected:
   int region_idx_;
   std::string ttree_root_file_;
   std::string rntuple_root_file_;

   static const std::vector<std::string> regions_;

   static std::string EnvOr(const char *name, const char *fallback)
   {
      const char *value = std::getenv(name);
      return (value && *value) ? value : fallback;
   }

   // Region queries print a per-call summary that would bury the benchmark output.
   // The old restore step reopened /dev/tty, which fails outright when the
   // benchmark runs from CI or with its output piped anywhere.
   void suppress_output()
   {
      fflush(stdout);
      saved_stdout_ = dup(fileno(stdout));
      null_ = fopen(NULL_DEVICE, "w");
      if (null_)
         dup2(fileno(null_), fileno(stdout));
   }

   void restore_output()
   {
      fflush(stdout);
      if (saved_stdout_ >= 0) {
         dup2(saved_stdout_, fileno(stdout));
         close(saved_stdout_);
         saved_stdout_ = -1;
      }
      if (null_) {
         fclose(null_);
         null_ = nullptr;
      }
   }

   int saved_stdout_ = -1;
   FILE *null_ = nullptr;

   const char *get_current_region() const { return regions_[region_idx_ % regions_.size()].c_str(); }
};

// Overridable so the set of regions matches whatever file the benchmark was given.
static std::vector<std::string> LoadRegions()
{
   const char *env = std::getenv("RAMTOOLS_BENCH_REGIONS");
   if (!env || !*env)
      return {"chr1:1000000-1001000", "chr1:1000000-2000000", "chr1:1-50000000", "chr21:1-48129895"};

   std::vector<std::string> regions;
   std::string spec = env;
   size_t start = 0;
   while (start <= spec.size()) {
      const size_t comma = spec.find(',', start);
      const std::string one = spec.substr(start, comma == std::string::npos ? std::string::npos : comma - start);
      if (!one.empty())
         regions.push_back(one);
      if (comma == std::string::npos)
         break;
      start = comma + 1;
   }
   return regions;
}

const std::vector<std::string> RegionQueryFixture::regions_ = LoadRegions();

BENCHMARK_DEFINE_F(RegionQueryFixture, TTree)(benchmark::State &state)
{
   if (ttree_root_file_.empty()) {
      state.SkipWithError("set RAMTOOLS_BENCH_TTREE to a file written by samtoram");
      return;
   }

   const char *region = get_current_region();
   int64_t total_reads_processed = 0;
   Long64_t reads_in_this_run = 0;

   for (auto _ : state) {
      suppress_output();
      reads_in_this_run = ramview(ttree_root_file_.c_str(), region, true, false, "perf.root");
      restore_output();

      total_reads_processed += reads_in_this_run;
   }

   state.SetItemsProcessed(total_reads_processed);
   state.counters["region_idx"] = region_idx_;
   state.SetLabel(std::to_string(reads_in_this_run) + " reads");
}

BENCHMARK_DEFINE_F(RegionQueryFixture, RNTuple)(benchmark::State &state)
{
   if (rntuple_root_file_.empty()) {
      state.SkipWithError("set RAMTOOLS_BENCH_RNTUPLE to a file written by samtoramntuple");
      return;
   }

   const char *region = get_current_region();
   int64_t total_reads_processed = 0;
   Long64_t reads_in_this_run = 0;

   for (auto _ : state) {
      suppress_output();
      reads_in_this_run = ramntupleview(rntuple_root_file_.c_str(), region, {true, false, "perf.root"});
      restore_output();

      total_reads_processed += reads_in_this_run;
   }

   state.SetItemsProcessed(total_reads_processed);
   state.counters["region_idx"] = region_idx_;
   state.SetLabel(std::to_string(reads_in_this_run) + " reads");
}

// One case per region. The indices used to be 0/3/6/9 into a fixed 18-region list;
// with a configurable list they wrapped and measured the same region twice.
BENCHMARK_REGISTER_F(RegionQueryFixture, TTree)->DenseRange(0, 3)->Unit(benchmark::kSecond);

BENCHMARK_REGISTER_F(RegionQueryFixture, RNTuple)->DenseRange(0, 3)->Unit(benchmark::kSecond);

BENCHMARK_MAIN();
