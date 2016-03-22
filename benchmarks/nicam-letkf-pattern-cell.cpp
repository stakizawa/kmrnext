/// \file
/// Benchmark program for NICAM-LETKF data access pattern
///
/// It requires KMR backend.
/// The number of MPI processes running this program should be
/// "kNumEnsemble x kNumRegion" because of the limitation of
/// DataStore::collate().
#include <iostream>
#include <sstream>
#include <cassert>
#include <ctime>
#include "kmrnext.hpp"

using namespace std;
using namespace kmrnext;

int rank = 0;

const size_t kNumIteration = 10;

const size_t kDimEnsembleData = 3;
//const size_t kDimRegionData   = 2;
const size_t kDimCellData     = 1;

#if DEBUG
const size_t kNumEnsemble  = 2;
const size_t kNumRegion    = 10;
const size_t kNumCell      = 10;
const size_t kElementCount = 2;
#else
const size_t kNumEnsemble  = 64;
const size_t kNumRegion    = 40;
const size_t kNumCell      = 1156;
// const size_t kNumCell      = 4624;
// const size_t kNumCell      = 18496;
// const size_t kNumCell      = 73984;

// Assume that each lattice has 6160 bytes of data (6160 = 1540 * 4)
const size_t kElementCount      = 1540;
#endif

const size_t kEnsembleDataDimSizes[kDimEnsembleData] =
  {kNumEnsemble, kNumRegion, kNumCell};

const bool kPrint = true;

class DataPrinter : public DataPack::Dumper {
public:
  string operator()(DataPack& dp)
  {
    ostringstream os;
    os << dp.key().to_string() << " : ";
    int *data = (int*)dp.data()->value();
    for (size_t i = 0; i < kElementCount; i++) {
      os << data[i] << " ";
    }
    os << endl;
    return os.str();
  }
};

struct Time {
  double loop_start;
  double loop_finish;

  double nicam_start;
  double nicam_finish;
  double letkf_start;
  double letkf_finish;

  double nicam_invoke;
  double nicam_cleanup;
  double letkf_invoke;
  double letkf_cleanup;

  double del0_start;
  double del0_finish;
  double del1_start;
  double del1_finish;

  double collate_start;
  double collate_finish;

  double loop() {
    return (loop_finish - loop_start) / 10E9;
  }
  double nicam() {
    return (nicam_finish - nicam_start) / 10E9;
  }
  double letkf() {
    return (letkf_finish - letkf_start) / 10E9;
  }
  double nicam_launch() {
    return (nicam_cleanup - nicam_invoke) / 10E9;
  }
  double letkf_launch() {
    return (letkf_cleanup - letkf_invoke) / 10E9;
  }
  double del0() {
    return (del0_finish - del0_start) / 10E9;
  }
  double del1() {
    return (del1_finish - del1_start) / 10E9;
  }
  double collate() {
    return (collate_finish - collate_start) / 10E9;
  }
};

void load_data(DataStore* ds);
void run_nicam(DataStore* inds, DataStore* outds, Time& time);
void run_letkf(DataStore* inds, DataStore* outds, Time& time);
void collate_ds(DataStore* ds, Time& time);
void print_line(string& str);
void print_line(ostringstream& os);
double gettime();
double gettime(DataStore::MapEnvironment& env);

//////////////////////////////////////////////////////////////////////////////
// Main starts from here.
//////////////////////////////////////////////////////////////////////////////
int
main(int argc, char **argv)
{
  KMRNext *next = KMRNext::init(argc, argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  {
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    if (nprocs != kNumEnsemble * kNumRegion) {
      if (rank == 0) {
	cout << "The number of processes should be "
	     << (kNumEnsemble * kNumRegion) << "." << endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  DataStore* ds0 = next->create_ds(kDimEnsembleData);
  ds0->set(kEnsembleDataDimSizes);
  load_data(ds0);
  DataPrinter dp;

  for (size_t i = 0; i < kNumIteration; i++) {
    Time time;
    time.loop_start = gettime();

    // run pseudo-NICAM
    DataStore* ds1 = next->create_ds(kDimEnsembleData);
    ds1->set(kEnsembleDataDimSizes);
    run_nicam(ds0, ds1, time);
    time.del0_start = gettime();
    delete ds0;
    time.del0_finish = gettime();

    // run pseudo-LETKF
    ds0 = next->create_ds(kDimEnsembleData);
    ds0->set(kEnsembleDataDimSizes);
    run_letkf(ds1, ds0, time);
    time.del1_start = gettime();
    delete ds1;
    time.del1_finish = gettime();

    // collate for parallel execution
    collate_ds(ds0, time);

    time.loop_finish = gettime();
    ostringstream os1;
    os1 << "Iteration[" << i << "]," << time.loop() << endl;
    os1 << "Invoke NICAM," << time.nicam_launch() << endl;
    os1 << "NICAM,"        << time.nicam() << endl;
    os1 << "Del NICAM In," << time.del0() << endl;
    os1 << "Invoke LETKF," << time.letkf_launch() << endl;
    os1 << "LETKF,"        << time.letkf() << endl;
    os1 << "Del LETKF In," << time.del1() << endl;
    os1 << "Collate,"      << time.collate() << endl;
    print_line(os1);
  }
  delete ds0;

  KMRNext::finalize();
  return 0;
}


class DataLoader : public DataStore::Loader<int> {
public:
  int operator()(DataStore* ds, const int& num)
  {
    Key key(kDimCellData);
    int *data_val = new int[kElementCount];
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (size_t i = 0; i < kElementCount; i++) {
      data_val[i] = num;
    }
    Data data((void*)data_val, sizeof(int) * kElementCount);

#ifdef _OPENMP
    #pragma omp parallel for firstprivate(key)
#endif
    for (size_t i = 0; i < kNumCell; i++) {
      key.set_dim(0, i);
      ds->add(key, data);
    }
    delete[] data_val;
    return 0;
  }
};

void load_data(DataStore* ds)
{
  vector<int> data_srcs;
  for (size_t i = 0; i < kNumEnsemble; i++) {
    for (size_t j = 0; j < kNumRegion; j++) {
      data_srcs.push_back((int)(j+1));
    }
  }
  DataLoader dl;
  ds->load_array(data_srcs, dl);
}

class PseudoNICAM : public DataStore::Mapper {
  Time& time_;
public:
  PseudoNICAM(Time& time) : time_(time) {};

  int operator()(DataStore* inds, DataStore* outds,
		 Key& key, vector<DataPack>& dps,
		 DataStore::MapEnvironment& env)
  {
#if DEBUG
    size_t local_count = dps.size();
    size_t total_count;
    MPI_Allreduce(&local_count, &total_count, 1, MPI_LONG, MPI_SUM,
		  env.mpi_comm);
    assert(total_count == (size_t)(kNumRegion * kNumCell));
#endif

    time_.nicam_start = gettime(env);
    for (vector<DataPack>::iterator itr = dps.begin(); itr != dps.end();
	 itr++) {
      int *data_new = new int[kElementCount];
      int *data_old = (int*)itr->data()->value();
#ifdef _OPENMP
      #pragma omp parallel for
#endif
      for (size_t i = 0; i < kElementCount; i++) {
	data_new[i] = data_old[i] + 1;
      }
      Data data(data_new, itr->data()->size());
      outds->add(itr->key(), data);
      delete[] data_new;
    }
    time_.nicam_finish = gettime(env);

    return 0;
  }
};

void run_nicam(DataStore* inds, DataStore* outds, Time& time)
{
  PseudoNICAM mapper(time);
  time.nicam_invoke = gettime();
  View view(kDimEnsembleData);
  bool view_flag[3] = {true, false, false};
  view.set(view_flag);
  inds->map(outds, mapper, view);
  time.nicam_cleanup = gettime();
}

class PseudoLETKF : public DataStore::Mapper {
  Time& time_;
public:
  PseudoLETKF(Time& time) : time_(time) {};

  int operator()(DataStore* inds, DataStore* outds,
		 Key& key, vector<DataPack>& dps,
		 DataStore::MapEnvironment& env)
  {
#if DEBUG
    int nprocs;
    MPI_Comm_size(env.mpi_comm, &nprocs);
    assert(nprocs == 1);
    size_t local_count = dps.size();
    assert(local_count == (size_t)kNumEnsemble);
#endif

    time_.letkf_start = gettime(env);
    for (vector<DataPack>::iterator itr = dps.begin(); itr != dps.end();
	 itr++) {
      int *data_new = new int[kElementCount];
      int *data_old = (int*)itr->data()->value();
#ifdef _OPENMP
      #pragma omp parallel for
#endif
      for (size_t i = 0; i < kElementCount; i++) {
	data_new[i] = data_old[i] - 1;
      }
      Data data(data_new, itr->data()->size());
      outds->add(itr->key(), data);
      delete[] data_new;
    }
    time_.letkf_finish = gettime(env);

    return 0;
  }
};

void run_letkf(DataStore* inds, DataStore* outds, Time& time)
{
  PseudoLETKF mapper(time);
  time.letkf_invoke = gettime();
  View view(kDimEnsembleData);
  bool view_flag[3] = {false, true, true};
  view.set(view_flag);
  inds->map_single(outds, mapper, view);
  time.letkf_cleanup = gettime();
}

void collate_ds(DataStore* ds, Time& time)
{
  time.collate_start = gettime();
  View view(kDimEnsembleData);
  bool view_flag[3] = {true, true, false};
  view.set(view_flag);
  ds->collate(view);
  time.collate_finish = gettime();
}

void print_line(string& str) {
  if (kPrint && (rank == 0)) {
    cout << str << endl;
  }
}

void print_line(ostringstream& os) {
  string s = os.str();
  print_line(s);
}

double gettime() {
  MPI_Barrier(MPI_COMM_WORLD);
  struct timespec ts;
  clock_gettime(CLOCK_REALTIME, &ts);
  return ((double) ts.tv_sec) * 10E9 + ((double) ts.tv_nsec);
}

double gettime(DataStore::MapEnvironment& env) {
  MPI_Barrier(env.mpi_comm);
  struct timespec ts;
  clock_gettime(CLOCK_REALTIME, &ts);
  return ((double) ts.tv_sec) * 10E9 + ((double) ts.tv_nsec);
}
