/// \file
/// Benchmark program for evaluating DataStore::map() operation
///
/// It requires KMR backend.
/// It measures time of DataStore::map() by changing the number of
/// nodes.  It is useful for evaluating scalability.  The allowed number
/// of nodes are 100, 200, 400, 800, 1600, 3200 and 6400.
///
/// The dimension of used DataStore is two, and its size dicided depending
/// on the number of nodes as follows.
///    - 100:  (100, 100)
///    - 200:  (200, 100)
///    - 400:  (400, 100)
///    - 800:  (800, 100)
///    - 1600: (1600, 100)
///    - 3200: (3200, 100)
///    - 6400: (6400, 100)
/// Used views are <F,F>, <F,T>, <T,F> and <T,T>.
///
/// In the DEBUG mode, the number of allowed processes are 2 and 4.  Sizes
/// of DataStores are (2, 100) and (4, 100).
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cassert>
#include <ctime>
#include "kmrnext.hpp"

using namespace std;
using namespace kmrnext;

#if DEBUG
const size_t kAllowedProcs   = 2;
const int    kNumProcList[2] = {2, 4};
const int    kNumIterations  = 2;
#else
const size_t kAllowedProcs   = 7;
const int    kNumProcList[7] = {100, 200, 400, 800, 1600, 3200, 6400};
const int    kNumIterations  = 10;
#endif

void setup(int argc, int* rankp, int* nprocsp);
void evaluate(KMRNext* next, int rank, int nprocs);


//////////////////////////////////////////////////////////////////////////////
// Main starts from here.
//////////////////////////////////////////////////////////////////////////////
int
main(int argc, char** argv)
{
  KMRNext *next = KMRNext::init(argc, argv);
  //next->enable_profile();
  int rank, nprocs;
  setup(argc, &rank, &nprocs);
  evaluate(next, rank, nprocs);
  KMRNext::finalize();
  return 0;
}


class Timer {
  vector<double> starts;
  vector<double> finishs;

  double gettime() {
    MPI_Barrier(MPI_COMM_WORLD);
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return ((double) ts.tv_sec) * 10E9 + ((double) ts.tv_nsec);
  }

public:
  void start() {
    starts.push_back(gettime());
  }

  void finish() {
    finishs.push_back(gettime());
  }

  string str(bool header=true, bool header_only=false) {
    assert(starts.size() == finishs.size());
    ostringstream os;
    if (header || header_only) {
      for (size_t i = 0; i < starts.size(); i++) {
	os << i << ",";
      }
      os << "average" << endl;
      if (header_only) {
	return os.str();
      }
    }
    double time_sum = 0;
    for (size_t i = 0; i < starts.size(); i++) {
      double time = finishs.at(i) - starts.at(i);
      time_sum += time;
      os << time << ",";
    }
    os << (time_sum / (double)starts.size()) << endl;
    return os.str();
  }
};

class DataLoader : public DataStore::Loader<long> {
  size_t element_count_;

public:
  DataLoader(size_t count) : element_count_(count) {}

  int operator()(DataStore* ds, const long& num)
  {
    size_t data_count = 1;
    for (size_t i = 0; i < ds->size(); i++) {
      data_count *= ds->dim(i);
    }

    long *value = new long[element_count_];
    for (size_t i = 0; i < element_count_; i++) {
      value[i] = num;
    }
    Data data((void*)value, sizeof(long) * element_count_);
    delete[] value;

    for (size_t i = 0; i < data_count; i++) {
      Key key = ds->index_to_key(i);
      ds->add(key, data);
    }
    return 0;
  }
};

class PseudoMapper : public DataStore::Mapper {
  int world_rank_;

public:
  PseudoMapper(int rank) : world_rank_(rank) {}

  int operator()(DataStore* inds, DataStore* outds,
		 Key& key, vector<DataPack>& dps,
		 DataStore::MapEnvironment& env)
  {
#if 0
#if DEBUG
    size_t local_count = dps.size();
    size_t total_count;
    MPI_Allreduce(&local_count, &total_count, 1, MPI_LONG, MPI_SUM,
		  env.mpi_comm);
    if (world_rank_ == 0) {
      cerr << "Total data count is " << total_count << endl;
    }
#endif
#endif
    return 0;
  }
};


void
setup(int argc, int* rankp, int* nprocsp)
{
  MPI_Comm_rank(MPI_COMM_WORLD, rankp);
  MPI_Comm_size(MPI_COMM_WORLD, nprocsp);
  bool included = false;
  for (size_t i = 0; i < kAllowedProcs; i++) {
    if (*nprocsp == kNumProcList[i]) {
      included = true;
      break;
    }
  }
  if (!included) {
    if (*rankp == 0) {
      cerr << "The Number of processes is illegal: " << (*nprocsp) << endl;
    }
    KMRNext::finalize();
    exit(1);
  }
}

void
evaluate(KMRNext* next, int rank, int nprocs)
{
  vector<long> datalist;
  for (size_t i = 0; i < (size_t)nprocs; i++) {
    datalist.push_back((long)(i+1));
  }
  DataLoader loader(1);

  size_t dim_ary[2] = {(size_t)nprocs, 100};
  DataStore *ds0 = next->create_ds(2);
  ds0->set(dim_ary);
  ds0->load_array(datalist, loader);

  bool ary_ff[2] = {false, false};
  View vff(2);
  vff.set(ary_ff);
  bool ary_ft[2] = {false, true};
  View vft(2);
  vft.set(ary_ft);
  bool ary_tf[2] = {true,  false};
  View vtf(2);
  vtf.set(ary_tf);
  bool ary_tt[2] = {true,  true};
  View vtt(2);
  vtf.set(ary_tt);

  Timer timer_ff;
  DataStore *ds = ds0->duplicate();
  for (size_t i = 0; i < kNumIterations; i++) {
    PseudoMapper mapper(rank);
    timer_ff.start();
    ds->map(mapper, vff);
    timer_ff.finish();
  }
  delete ds;

  Timer timer_ft;
  ds = ds0->duplicate();
  for (size_t i = 0; i < kNumIterations; i++) {
    PseudoMapper mapper(rank);
    timer_ft.start();
    ds->map(mapper, vft);
    timer_ft.finish();
  }
  delete ds;

  Timer timer_tf;
  ds = ds0->duplicate();
  for (size_t i = 0; i < kNumIterations; i++) {
    PseudoMapper mapper(rank);
    timer_tf.start();
    ds->map(mapper, vtf);
    timer_tf.finish();
  }
  delete ds;

  Timer timer_tt;
  ds = ds0->duplicate();;
  for (size_t i = 0; i < kNumIterations; i++) {
    PseudoMapper mapper(rank);
    timer_tt.start();
    ds->map(mapper, vtt);
    timer_tt.finish();
  }
  delete ds;
  delete ds0;
  if (rank == 0) {
    cout << "  ," << timer_ff.str(true, true);
    cout << "FF," << timer_ff.str(false);
    cout << "FT," << timer_ft.str(false);
    cout << "TF," << timer_tf.str(false);
    cout << "TT," << timer_tt.str(false) << endl;
  }
}
