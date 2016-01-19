/// \file
/// Benchmark program for NICAM-LETKF data access pattern
///
/// It requires KMR.
/// It is recomended that the number of MPI processes running this program
/// is kNumEnsemble x kNumRegion as this configuration can achieve a good
/// performance with the least number of processes.
#include <iostream>
#include <sstream>
#include <cassert>
#include <ctime>
#include <unistd.h>
#include "kmrnext.hpp"

using namespace std;
using namespace kmrnext;

int rank = 0;

const size_t kNumIteration = 20;

const size_t kDimEnsembleData = 3;
//const size_t kDimRegionData   = 2;
const size_t kDimCellData     = 1;

#if DEBUG
const size_t kNumEnsemble	= 2;
const size_t kNumRegion		= 10;
const size_t kNumCell		= 10;
const size_t kElementCount	= 2;
const unsigned int kTimeNICAM	= 0; // msec
const unsigned int kTimeLETKF	= 0; // msec
#else
const size_t kNumEnsemble       = 64;
const size_t kNumRegion         = 10;
const size_t kNumCell           = 1156;
// Assume that each lattice has 6160 KB of data (6160 = 1540 * 4)
const size_t kElementCount      = 1540;
const unsigned int kTimeNICAM   = 50000; // msec
const unsigned int kTimeLETKF   = 5000;  // msec
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
};

void load_data(DataStore* ds);
void run_nicam(DataStore* inds, DataStore* outds, Time& time);
void run_letkf(DataStore* inds, DataStore* outds, Time& time);
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

  DataStore* ds0 = next->create_ds(kDimEnsembleData);
  ds0->set(kEnsembleDataDimSizes);
  load_data(ds0);
  DataPrinter dp;

  for (size_t i = 0; i < kNumIteration; i++) {
    ostringstream os0;
    os0 << "Iteration[" << i << "] starts.";
    print_line(os0);
    Time time;
    time.loop_start = gettime();

    // run pseudo-NICAM
    DataStore* ds1 = next->create_ds(kDimEnsembleData);
    ds1->set(kEnsembleDataDimSizes);
    run_nicam(ds0, ds1, time);
    delete ds0;

    // run pseudo-LETKF
    ds0 = next->create_ds(kDimEnsembleData);
    ds0->set(kEnsembleDataDimSizes);
    run_letkf(ds1, ds0, time);
    delete ds1;

    time.loop_finish = gettime();
    ostringstream os1;
    os1 << "Iteration[" << i << "] ends in " << time.loop() << " sec." << endl;
    os1 << "  Invoking NICAM takes " << time.nicam_launch() << " sec." << endl;
    os1 << "    NICAM consumes " << time.nicam() << " sec." << endl;
    os1 << "  Invoking LETKF takes " << time.letkf_launch() << " sec." << endl;
    os1 << "    LETKF consumes " << time.letkf() << " sec." << endl;
    print_line(os1);
  }

  KMRNext::finalize();
  return 0;
}


class DataLoader : public DataStore::Loader<int> {
public:
  int operator()(DataStore* ds, const int& num)
  {
    Key key(kDimCellData);
    int *data_val = new int[kElementCount];
    for (size_t i = 0; i < kElementCount; i++) {
      data_val[i] = num;
    }
    Data data((void*)data_val, sizeof(int) * kElementCount);

    for (size_t i = 0; i < kNumCell; i++) {
      key.set_dim(0, i);
      ds->add(key, data);
    }
    delete data_val;
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
    usleep(kTimeNICAM);
    time_.nicam_finish = gettime(env);

    for (vector<DataPack>::iterator itr = dps.begin(); itr != dps.end();
	 itr++) {
      int *data_new = new int[kElementCount];
      int *data_old = (int*)itr->data()->value();
      for (size_t i = 0; i < kElementCount; i++) {
	data_new[i] = data_old[i] + 1;
      }
      Data data(data_new, itr->data()->size());
      outds->add(itr->key(), data);
      delete data_new;
    }
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
    size_t local_count = dps.size();
    size_t total_count;
    MPI_Allreduce(&local_count, &total_count, 1, MPI_LONG, MPI_SUM,
		  env.mpi_comm);
    assert(total_count == (size_t)kNumEnsemble);
#endif

    time_.letkf_start = gettime(env);
    int  sndcnt = (int)(dps.size() * kElementCount);
    int *sndbuf = new int[sndcnt];
    int *rcvbuf = new int[kNumEnsemble * kElementCount];
    size_t sndbuf_idx = 0;
    for (vector<DataPack>::iterator itr = dps.begin(); itr != dps.end();
	 itr++) {
      int *data = (int*)itr->data()->value();
      size_t dsize = itr->data()->size();
      memcpy(sndbuf + sndbuf_idx, data, dsize);
      sndbuf_idx += dsize;
    }
    MPI_Allgather(sndbuf, sndcnt, MPI_INT, rcvbuf, sndcnt, MPI_INT,
		  env.mpi_comm);
    delete sndbuf;
    delete rcvbuf;
    usleep(kTimeLETKF);
    time_.letkf_finish = gettime(env);

    for (vector<DataPack>::iterator itr = dps.begin(); itr != dps.end();
	 itr++) {
      int *data_new = new int[kElementCount];
      int *data_old = (int*)itr->data()->value();
      for (size_t i = 0; i < kElementCount; i++) {
	data_new[i] = data_old[i] - 1;
      }
      Data data(data_new, itr->data()->size());
      outds->add(itr->key(), data);
      delete data_new;
    }
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
  inds->map(outds, mapper, view);
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
  return ((double) ts.tv_sec) * 10E9 + ((double) ts.tv_nsec);
}

double gettime(DataStore::MapEnvironment& env) {
  MPI_Barrier(env.mpi_comm);
  struct timespec ts;
  clock_gettime(CLOCK_REALTIME, &ts);
  return ((double) ts.tv_sec) * 10E9 + ((double) ts.tv_nsec);
}
