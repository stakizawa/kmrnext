/// \file
/// Benchmark program for evaluating DataStore::map() operation
///
/// It requires KMR backend.
/// It measures time for DataStore::map() by changing the following
/// parameters; # of dimensions, View type, # of data count in a DataStore
/// and size of each data.  This program accepts the following options to
/// perform the evaluation.
///    - dimension
///    - view
///    - count
///    - size
///
/// The number of processes to run this program should be two for DEBUG mode
/// and 1000 for large scale run.
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cassert>
#include <ctime>
#include "kmrnext.hpp"

using namespace std;
using namespace kmrnext;

#if DEBUG
const int kNumProcs      = 2;
const int kNumIterations = 2;
#else
const int kNumProcs      = 1000;
const int kNumIterations = 10;
#endif

void setup(int argc, int* rankp, int* nprocsp);
void eval_dimension(KMRNext* next, int rank, int nprocs);
void eval_view(KMRNext* next, int rank, int nprocs);
void eval_count(KMRNext* next, int rank, int nprocs);
void eval_size(KMRNext* next, int rank, int nprocs);
void help();


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
  string option(argv[1]);
  if (option == "dimension") {
    eval_dimension(next, rank, nprocs);
  } else if (option == "view") {
    eval_view(next, rank, nprocs);
  } else if (option == "count") {
    eval_count(next, rank, nprocs);
  } else if (option == "size") {
    eval_size(next, rank, nprocs);
  } else {
    help();
    return 0;
  }

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

  string str(bool header=true) {
    assert(starts.size() == finishs.size());
    ostringstream os;
    if (header) {
      for (size_t i = 0; i < starts.size(); i++) {
	os << i << ",";
      }
      os << endl;
    }
    for (size_t i = 0; i < starts.size(); i++) {
      double time = finishs.at(i) - starts.at(i);
      os << time << ",";
    }
    os << endl;
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
#if DEBUG
    size_t local_count = dps.size();
    size_t total_count;
    MPI_Allreduce(&local_count, &total_count, 1, MPI_LONG, MPI_SUM,
		  env.mpi_comm);
    if (world_rank_ == 0) {
      cerr << "Total data count is " << total_count << endl;
    }
#endif
    return 0;
  }
};


void
setup(int argc, int* rankp, int* nprocsp)
{
  if (argc != 2) {
    help();
    exit(1);
  }
  MPI_Comm_rank(MPI_COMM_WORLD, rankp);
  MPI_Comm_size(MPI_COMM_WORLD, nprocsp);
  if (*nprocsp < kNumProcs) {
    if (*rankp == 0) {
      cerr << "Number of processes should be bigger than " << kNumProcs
	   << endl;
    }
    KMRNext::finalize();
    exit(1);
  }
  if (*nprocsp != kNumProcs) {
    if (*rankp == 0) {
      cerr << "Number of runnning processes is bigger than the required "
	   << "number of processes." << endl;
      cerr << (*nprocsp - kNumProcs) << " processes are not used." << endl;
    }
  }
}

static void
eval_dimension_once(KMRNext* next, int rank, int ndim, size_t* dim_ary,
		    bool* view_ary, vector<long>& datalist, DataLoader& loader)
{
  DataStore *ds = next->create_ds(ndim);
  ds->set(dim_ary);
  ds->load_array(datalist, loader);
  View view(ndim);
  view.set(view_ary);
  Timer timer;
  for (size_t i = 0; i < kNumIterations; i++) {
    PseudoMapper mapper(rank);
    timer.start();
    ds->map(NULL, mapper, view);
    timer.finish();
  }
  delete ds;
  cout << timer.str() << endl;
}

void
eval_dimension(KMRNext* next, int rank, int nprocs)
{
#if DEBUG
  size_t ary2[2] = {2, 10000};
  size_t ary3[3] = {2, 10, 1000};
  size_t ary4[4] = {2, 10, 10, 100};
  size_t ary5[5] = {2, 10, 10, 10, 10};
  size_t ary6[6] = {2, 10, 10, 10, 5, 2};
  size_t ary7[7] = {2, 10, 10, 5, 2, 5, 2};
  size_t ary8[8] = {2, 10, 5, 2, 5, 2, 5, 2};
#else
  size_t ary2[2] = {1000, 10000};
  size_t ary3[3] = {1000, 10, 1000};
  size_t ary4[4] = {1000, 10, 10, 100};
  size_t ary5[5] = {1000, 10, 10, 10, 10};
  size_t ary6[6] = {1000, 10, 10, 10, 5, 2};
  size_t ary7[7] = {1000, 10, 10, 5, 2, 5, 2};
  size_t ary8[8] = {1000, 10, 5, 2, 5, 2, 5, 2};
#endif

  bool   aryv2[2] = {true, false};
  bool   aryv3[3] = {true, false, false};
  bool   aryv4[4] = {true, false, false, false};
  bool   aryv5[5] = {true, false, false, false, false};
  bool   aryv6[6] = {true, false, false, false, false, false};
  bool   aryv7[7] = {true, false, false, false, false, false, false};
  bool   aryv8[8] = {true, false, false, false, false, false, false, false};

  vector<long> datalist;
  for (size_t i = 0; i < kNumProcs; i++) {
    datalist.push_back((long)(i+1));
  }
  DataLoader loader(1);

  eval_dimension_once(next, rank, 2, ary2, aryv2, datalist, loader);
  eval_dimension_once(next, rank, 3, ary3, aryv3, datalist, loader);
  eval_dimension_once(next, rank, 4, ary4, aryv4, datalist, loader);
  eval_dimension_once(next, rank, 5, ary5, aryv5, datalist, loader);
  eval_dimension_once(next, rank, 6, ary6, aryv6, datalist, loader);
  eval_dimension_once(next, rank, 7, ary7, aryv7, datalist, loader);
  eval_dimension_once(next, rank, 8, ary8, aryv8, datalist, loader);
}

void
eval_view(KMRNext* next, int rank, int nprocs)
{
  throw runtime_error("not implemented yet.");
}

void
eval_count(KMRNext* next, int rank, int nprocs)
{
  throw runtime_error("not implemented yet.");
}

void
eval_size(KMRNext* next, int rank, int nprocs)
{
  throw runtime_error("not implemented yet.");
}

void
help()
{
  cout << "Options: " << endl;
  cout << "  dimension" << endl;
  cout << "      Change number of dimensions of a DataStore" << endl;
  cout << "  view" << endl;
  cout << "      Change view type" << endl;
  cout << "  count" << endl;
  cout << "      Change count of data in a DataStore" << endl;
  cout << "  size" << endl;
  cout << "      Change size of data in a DataStore" << endl;
  cout << endl;
  KMRNext::finalize();
}
