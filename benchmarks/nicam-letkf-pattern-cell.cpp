/// \file
/// Benchmark program for NICAM-LETKF data access pattern
///
/// It requires KMR backend.
/// The number of MPI processes running this program should be
/// "kNumEnsemble x kNumRegion".
#include <iostream>
#include <sstream>
#include <cassert>
#include <ctime>
#include <unistd.h>
#include "kmrnext.hpp"

using namespace std;
using namespace kmrnext;

int rank = 0;
int nprocs = 1;

// If true, map works in-place.
const bool kMapInplace = true;

// If true, set the best split
const bool kSetSplit = true;

// If true, output in YAML format
const bool kOutputYAML = true;

const size_t kNumIteration = 10;

const size_t kDimEnsembleData = 3;

#if DEBUG
const size_t kNumEnsemble     = 2;
const size_t kNumRegion       = 10;
const size_t kNumCell         = 10;
const size_t kElementCount    = 2;
const unsigned int kTimeNICAM = 0; // msec
const unsigned int kTimeLETKF = 0; // msec
#else
const size_t kNumEnsemble     = 64;
const size_t kNumRegion       = 40;
const size_t kNumCell         = 1156;
// const size_t kNumCell         = 4624;
// const size_t kNumCell         = 18496;
// const size_t kNumCell         = 73984;

// Assume that each lattice has 6160 bytes of data (6160 = 1540 * 4)
const size_t kElementCount    = 1540;
const unsigned int kTimeNICAM = 94000; // msec
const unsigned int kTimeLETKF = 1000;   // msec
#endif

const size_t kEnsembleDataDimSizes[kDimEnsembleData] =
  {kNumEnsemble, kNumRegion, kNumCell};

const bool kPrint = true;

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

  double alc0_start;
  double alc0_finish;
  double alc1_start;
  double alc1_finish;

  double del0_start;
  double del0_finish;
  double del1_start;
  double del1_finish;

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
  double alc0() {
    return (alc0_finish - alc0_start) / 10E9;
  }
  double alc1() {
    return (alc1_finish - alc1_start) / 10E9;
  }
  double del0() {
    return (del0_finish - del0_start) / 10E9;
  }
  double del1() {
    return (del1_finish - del1_start) / 10E9;
  }
};

void load_data(DataStore* ds);
void run_nicam(DataStore* inds, DataStore* outds, Time& time);
void run_letkf(DataStore* inds, DataStore* outds, Time& time);
void print_line(string& str);
void print_line(ostringstream& os);
double gettime();
double gettime(DataStore::MapEnvironment& env);
int calculate_task_nprocs(View& view, View& alc_view, int given_nprocs);

//////////////////////////////////////////////////////////////////////////////
// Main starts from here.
//////////////////////////////////////////////////////////////////////////////
int
main(int argc, char **argv)
{
  KMRNext *next = KMRNext::init(argc, argv);
  next->set_io_mode(kmrnext::KMRNext::File);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  {
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

  double *loop_times = new double[kNumIteration];
  DataStore* ds1;
  for (size_t i = 0; i < kNumIteration; i++) {
    Time time;
    time.loop_start = gettime();

    // run pseudo-NICAM
    if (kMapInplace) {
      ds1 = ds0;
      time.alc0_start = gettime();
      time.alc0_finish = gettime();
      run_nicam(ds0, ds1, time);
      time.del0_start = gettime();
      time.del0_finish = gettime();
    } else {
      time.alc0_start = gettime();
      ds1 = next->create_ds(kDimEnsembleData);
      ds1->set(kEnsembleDataDimSizes);
      time.alc0_finish = gettime();
      run_nicam(ds0, ds1, time);
      time.del0_start = gettime();
      delete ds0;
      time.del0_finish = gettime();
    }

    // run pseudo-LETKF
    if (kMapInplace) {
      ds0 = ds1;
      time.alc1_start = gettime();
      time.alc1_finish = gettime();
      run_letkf(ds1, ds0, time);
      time.del1_start = gettime();
      time.del1_finish = gettime();
    } else {
      time.alc1_start = gettime();
      ds0 = next->create_ds(kDimEnsembleData);
      ds0->set(kEnsembleDataDimSizes);
      time.alc1_finish = gettime();
      run_letkf(ds1, ds0, time);
      time.del1_start = gettime();
      delete ds1;
      time.del1_finish = gettime();
    }

    time.loop_finish = gettime();
    loop_times[i] = time.loop();
    if (!kOutputYAML) {
      ostringstream os1;
      os1 << "Iteration[" << i << "]," << time.loop() << endl;
      os1 << "Alc NICAM In," << time.alc0() << endl;
      os1 << "Invoke NICAM," << time.nicam_launch() << endl;
      //os1 << "NICAM,"        << time.nicam() << endl;
      os1 << "Del NICAM In," << time.del0() << endl;
      os1 << "Alc LETKF In," << time.alc1() << endl;
      os1 << "Invoke LETKF," << time.letkf_launch() << endl;
      //os1 << "LETKF,"        << time.letkf() << endl;
      os1 << "Del LETKF In," << time.del1() << endl;
      print_line(os1);
    }
  }
  delete ds0;

  if (kOutputYAML) {
    ostringstream os0;
    for (size_t i = 0; i < kNumIteration; i++) {
      os0 << "  - " << loop_times[i] << endl;
    }
    print_line(os0);
  }
  delete[] loop_times;

  KMRNext::finalize();
  return 0;
}


class DataLoader : public DataStore::Loader<long> {
public:
  int operator()(DataStore* ds, const long& num)
  {
    size_t x = num / kNumRegion;
    size_t y = num % kNumRegion;
    int *data_val = new int[kElementCount];
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (size_t i = 0; i < kElementCount; i++) {
      data_val[i] = static_cast<int>(y) + 1;
    }
    Data data(static_cast<void*>(data_val), sizeof(int) * kElementCount);

    Key key(kDimEnsembleData);
    key.set_dim(0, x);
    key.set_dim(1, y);
#ifdef _OPENMP
    #pragma omp parallel for firstprivate(key)
#endif
    for (size_t i = 0; i < kNumCell; i++) {
      key.set_dim(2, i);
      ds->add(key, data);
    }
    delete[] data_val;
    return 0;
  }
};

void load_data(DataStore* ds)
{
  vector<long> data_srcs;
  for (size_t i = 0; i < kNumEnsemble; i++) {
    for (size_t j = 0; j < kNumRegion; j++) {
      data_srcs.push_back(i * kNumRegion + j);
    }
  }
  DataLoader dl;
  ds->load_integers(data_srcs, dl);
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
    int nprocs_nicam;
    MPI_Comm_size(env.mpi_comm, &nprocs_nicam);
    int nprocs_calc = calculate_task_nprocs(env.view, env.split,
					    nprocs_nicam);
    assert(nprocs_nicam == nprocs_calc);
    size_t local_count = dps.size();
    size_t total_count;
    MPI_Allreduce(&local_count, &total_count, 1, MPI_LONG, MPI_SUM,
		  env.mpi_comm);
    assert(total_count == static_cast<size_t>(kNumRegion * kNumCell));
#endif

    time_.nicam_start = gettime(env);
    usleep(kTimeNICAM * 1000);
    for (vector<DataPack>::iterator itr = dps.begin(); itr != dps.end();
	 itr++) {
      int *data_new = new int[kElementCount];
      int *data_old = static_cast<int*>(itr->data().value());
#ifdef _OPENMP
      #pragma omp parallel for
#endif
      for (size_t i = 0; i < kElementCount; i++) {
	data_new[i] = data_old[i] + 1;
      }
      Data data(data_new, itr->data().size());
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

  if (kSetSplit) {
    View split(kDimEnsembleData);
    bool split_flag[3] = {true, true, false};
    split.set(split_flag);
    inds->set_split(split);
  }

  View view(kDimEnsembleData);
  bool view_flag[3] = {true, false, false};
  view.set(view_flag);
  inds->map(mapper, view, outds);

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
    int nprocs_letkf;
    MPI_Comm_size(env.mpi_comm, &nprocs_letkf);
    assert(nprocs_letkf == 1);
    int nprocs_calc = calculate_task_nprocs(env.view, env.split,
					    nprocs_letkf);
    assert(nprocs_letkf == nprocs_calc);
    size_t local_count = dps.size();
    assert(local_count == static_cast<size_t>(kNumEnsemble));
#endif

    time_.letkf_start = gettime(env);
    usleep(kTimeLETKF * 1000);
    for (vector<DataPack>::iterator itr = dps.begin(); itr != dps.end();
	 itr++) {
      int *data_new = new int[kElementCount];
      int *data_old = static_cast<int*>(itr->data().value());
#ifdef _OPENMP
      #pragma omp parallel for
#endif
      for (size_t i = 0; i < kElementCount; i++) {
	data_new[i] = data_old[i] - 1;
      }
      Data data(data_new, itr->data().size());
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

  if (kSetSplit) {
    View split(kDimEnsembleData);
    bool split_flag[3] = {false, true, true};
    split.set(split_flag);
    inds->set_split(split);
  }

  View view(kDimEnsembleData);
  bool view_flag[3] = {false, true, true};
  view.set(view_flag);
  inds->map(mapper, view, outds);

  time.letkf_cleanup = gettime();
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
  return (static_cast<double>(ts.tv_sec) * 10E9 +
	  static_cast<double>(ts.tv_nsec));
}

double gettime(DataStore::MapEnvironment& env) {
  MPI_Barrier(env.mpi_comm);
  struct timespec ts;
  clock_gettime(CLOCK_REALTIME, &ts);
  return (static_cast<double>(ts.tv_sec) * 10E9 +
	  static_cast<double>(ts.tv_nsec));
}

int calculate_task_nprocs(View& view, View& alc_view, int given_nprocs) {
  int total_nprocs = 1;
  int nprocs_calc = 1;
  for (size_t i = 0; i < view.size(); i++) {
    if (!view.dim(i) && alc_view.dim(i)) {
      nprocs_calc *= static_cast<int>(kEnsembleDataDimSizes[i]);
    }
    if (alc_view.dim(i)) {
      total_nprocs *= static_cast<int>(kEnsembleDataDimSizes[i]);
    }
  }
  if (nprocs < total_nprocs) {
    nprocs_calc = given_nprocs;
  }
  return nprocs_calc;
}
