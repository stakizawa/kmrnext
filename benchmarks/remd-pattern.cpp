/// \file
/// Benchmark program for REMD data access pattern
///
/// It requires KMR backend.
/// It is recomended that the number of MPI processes running this program
/// is kNumEnsemble x kNumProc as this configuration can achieve a good
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

const size_t kNumIteration = 10;

const size_t kNumDimensions = 3;
#if DEBUG
const size_t kNumEnsemble   = 2;
const size_t kNumProc       = 8;
const size_t kNumData       = 6;
const size_t kElementCount  = 1;    // a double value
const unsigned int kTimeMD  = 1000; // msec
const unsigned int kTimeEX  = 1;    // msec
#else
const size_t kNumEnsemble   = 128;
const size_t kNumProc       = 8;
const size_t kNumData       = 6;
const size_t kElementCount  = 1;     // a double value
const unsigned int kTimeMD  = 20000; // msec
const unsigned int kTimeEX  = 1   ;  // msec
#endif

const size_t kDataStoreSizes[kNumDimensions] =
  {kNumEnsemble, kNumProc, kNumData};

const bool kPrint = true;

#if DEBUG
class DataPrinter : public DataPack::Dumper {
public:
  string operator()(DataPack& dp)
  {
    ostringstream os;
    os << dp.key().to_string() << " : ";
    double *data = (double*)dp.data()->value();
    for (size_t i = 0; i < kElementCount; i++) {
      os << data[i] << " ";
    }
    os << endl;
    return os.str();
  }
};
#endif

struct Time {
  double loop_start;
  double loop_finish;

  double md_start;
  double md_finish;
  double ex_start;
  double ex_finish;

  double md_invoke;
  double md_cleanup;
  double ex_invoke;
  double ex_cleanup;

  double loop() {
    return (loop_finish - loop_start) / 10E9;
  }
  double md() {
    return (md_finish - md_start) / 10E9;
  }
  double ex() {
    return (ex_finish - ex_start) / 10E9;
  }
  double md_launch() {
    return (md_cleanup - md_invoke) / 10E9;
  }
  double ex_launch() {
    return (ex_cleanup - ex_invoke) / 10E9;
  }
};

void load_data(DataStore* ds);
void run_md(DataStore* inds, DataStore* outds, Time& time);
void run_ex(DataStore* inds, DataStore* outds, Time& time);
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

  DataStore* ds0 = next->create_ds(kNumDimensions);
  ds0->set(kDataStoreSizes);
  load_data(ds0);

  for (size_t i = 0; i < kNumIteration; i++) {
    ostringstream os0;
    os0 << "Iteration[" << i << "] starts.";
    print_line(os0);
    Time time;
    time.loop_start = gettime();

    // run pseudo-MD
    DataStore* ds1 = next->create_ds(kNumDimensions);
    ds1->set(kDataStoreSizes);
    run_md(ds0, ds1, time);
    delete ds0;

    // run pseudo-EX
    ds0 = next->create_ds(kNumDimensions);
    ds0->set(kDataStoreSizes);
    run_ex(ds1, ds0, time);
    delete ds1;

    time.loop_finish = gettime();
    ostringstream os1;
    os1 << "Iteration[" << i << "] ends in " << time.loop() << " sec." << endl;
    os1 << "  Invoking MD takes " << time.md_launch() << " sec." << endl;
    os1 << "    MD consumes " << time.md() << " sec." << endl;
    os1 << "  Invoking EX takes " << time.ex_launch() << " sec." << endl;
    os1 << "    EX consumes " << time.ex() << " sec." << endl;
    print_line(os1);
  }

  KMRNext::finalize();
  return 0;
}


class DataLoader : public DataStore::Loader<double> {
public:
  int operator()(DataStore* ds, const double& num)
  {
    Key key(1);
    double *data_val = new double[kElementCount];
    for (size_t i = 0; i < kElementCount; i++) {
      data_val[i] = num;
    }
    Data data((void*)data_val, sizeof(double) * kElementCount);

    for (size_t i = 0; i < kNumData; i++) {
      key.set_dim(0, i);
      ds->add(key, data);
    }
    delete[] data_val;
    return 0;
  }
};

void load_data(DataStore* ds)
{
  vector<double> data_srcs;
  for (size_t i = 0; i < kNumEnsemble; i++) {
    for (size_t j = 0; j < kNumProc; j++) {
      data_srcs.push_back((double)(i+1));
    }
  }
  DataLoader dl;
  ds->load_array(data_srcs, dl);
}

class PseudoMD : public DataStore::Mapper {
  Time& time_;
public:
  PseudoMD(Time& time) : time_(time) {};

  int operator()(DataStore* inds, DataStore* outds,
		 Key& key, vector<DataPack>& dps,
		 DataStore::MapEnvironment& env)
  {
#if DEBUG
    size_t local_count = dps.size();
    size_t total_count;
    MPI_Allreduce(&local_count, &total_count, 1, MPI_LONG, MPI_SUM,
		  env.mpi_comm);
    assert(total_count == (size_t)(kNumProc * kNumData));
#endif

    time_.md_start = gettime(env);
    usleep(kTimeMD * 1000);
    time_.md_finish = gettime(env);

    for (vector<DataPack>::iterator itr = dps.begin(); itr != dps.end();
	 itr++) {
      double *data_new = new double[kElementCount];
      double *data_old = (double*)itr->data()->value();
      for (size_t i = 0; i < kElementCount; i++) {
	data_new[i] = data_old[i] + 1;
      }
      Data data(data_new, itr->data()->size());
      outds->add(itr->key(), data);
      delete[] data_new;
    }
    return 0;
  }
};

void run_md(DataStore* inds, DataStore* outds, Time& time)
{
  PseudoMD mapper(time);
  time.md_invoke = gettime();
  View view(kNumDimensions);
  bool view_flag[3] = {true, false, false};
  view.set(view_flag);
  inds->map(outds, mapper, view);
  time.md_cleanup = gettime();
}

class PseudoEX : public DataStore::Mapper {
  Time& time_;
public:
  PseudoEX(Time& time) : time_(time) {};

  int operator()(DataStore* inds, DataStore* outds,
		 Key& key, vector<DataPack>& dps,
		 DataStore::MapEnvironment& env)
  {
#if DEBUG
    size_t local_count = dps.size();
    size_t total_count;
    MPI_Allreduce(&local_count, &total_count, 1, MPI_LONG, MPI_SUM,
		  env.mpi_comm);
    assert(total_count == (size_t)(kNumEnsemble * kNumData));
#endif

    time_.ex_start = gettime(env);
    // perform allgather to collect all ensemble data
    int ex_nprocs;
    int ex_rank;
    MPI_Comm_size(env.mpi_comm, &ex_nprocs);
    MPI_Comm_rank(env.mpi_comm, &ex_rank);
    int nensemble = (int)(dps.size() / kNumData);
    int *nensembles = new int[ex_nprocs];
    MPI_Allgather(&nensemble, 1, MPI_INT, nensembles, 1, MPI_INT,
		  env.mpi_comm);
    // setup send data
    int     sndcnt = (int)(dps.size() * kElementCount);
    double *sndbuf = new double[sndcnt];
    size_t sndbuf_idx = 0;
    for (vector<DataPack>::iterator itr = dps.begin(); itr != dps.end();
	 itr++) {
      double *data = (double*)itr->data()->value();
      size_t dsize = itr->data()->size();
      memcpy(sndbuf + sndbuf_idx, data, dsize);
      sndbuf_idx += (dsize / sizeof(double));
    }
    // setup recv data
    int *rcvcnts = new int[ex_nprocs];
    int *rdispls = new int[ex_nprocs];
    int rcvbuf_siz = 0;
    for (int i = 0; i < ex_nprocs; i++) {
      rcvcnts[i] = nensembles[i] * (int)(kNumData * kElementCount);
      rcvbuf_siz += rcvcnts[i];
      if (i == 0) {
	rdispls[i] = 0;
      } else {
	rdispls[i] = rdispls[i-1] + rcvcnts[i-1];
      }
    }
    assert(rcvbuf_siz == (int)(kNumEnsemble * kNumData * kElementCount));
    double *rcvbuf = new double[rcvbuf_siz];
    // perform allgatherv
    MPI_Allgatherv(sndbuf, sndcnt, MPI_DOUBLE, rcvbuf, rcvcnts, rdispls,
		   MPI_DOUBLE, env.mpi_comm);
    delete[] rcvbuf;
    delete[] rdispls;
    delete[] rcvcnts;
    delete[] sndbuf;
    delete[] nensembles;
    usleep(kTimeEX * 1000);
    time_.ex_finish = gettime(env);

    for (vector<DataPack>::iterator itr = dps.begin(); itr != dps.end();
	 itr++) {
      double *data_new = new double[kElementCount];
      double *data_old = (double*)itr->data()->value();
      for (size_t i = 0; i < kElementCount; i++) {
	data_new[i] = data_old[i] - 1;
      }
      Data data(data_new, itr->data()->size());
      outds->add(itr->key(), data);
      delete[] data_new;
    }
    return 0;
  }
};

void run_ex(DataStore* inds, DataStore* outds, Time& time)
{
  PseudoEX mapper(time);
  time.ex_invoke = gettime();
  View view(kNumDimensions);
  bool view_flag[3] = {false, true, false};
  view.set(view_flag);
  inds->map(outds, mapper, view);
  time.ex_cleanup = gettime();
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
#ifdef BACKEND_KMR
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  struct timespec ts;
  clock_gettime(CLOCK_REALTIME, &ts);
  return ((double) ts.tv_sec) * 10E9 + ((double) ts.tv_nsec);
}

double gettime(DataStore::MapEnvironment& env) {
#ifdef BACKEND_KMR
  MPI_Barrier(env.mpi_comm);
#endif
  struct timespec ts;
  clock_gettime(CLOCK_REALTIME, &ts);
  return ((double) ts.tv_sec) * 10E9 + ((double) ts.tv_nsec);
}
