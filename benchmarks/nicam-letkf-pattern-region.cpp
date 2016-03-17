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
const size_t kNumRegion    = 10;
const size_t kNumCell      = 1156;
// const size_t kNumCell      = 4624;
// const size_t kNumCell      = 18496;
// const size_t kNumCell      = 73984;

// Assume that each grid point has 6160 bytes of data (6160 = 1540 * 4)
const size_t kElementCount = 1540;
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
#if 0
  string str0 = ds0->dump(dp);
  print_line(str0);
#endif

  for (size_t i = 0; i < kNumIteration; i++) {
    Time time;
    time.loop_start = gettime();

    // run pseudo-NICAM
    DataStore* ds1 = next->create_ds(kDimEnsembleData);
    ds1->set(kEnsembleDataDimSizes);
    run_nicam(ds0, ds1, time);
    delete ds0;
#if 0
    string str1 = ds1->dump(dp);
    print_line(str1);
#endif

    // run pseudo-LETKF
    ds0 = next->create_ds(kDimEnsembleData);
    ds0->set(kEnsembleDataDimSizes);
    run_letkf(ds1, ds0, time);
    delete ds1;
#if 0
    string str2 = ds0->dump(dp);
    print_line(str2);
#endif

    time.loop_finish = gettime();
    ostringstream os1;
    os1 << "Iteration[" << i << "]," << time.loop() << endl;
    os1 << "Invoke NICAM," << time.nicam_launch() << endl;
    os1 << "NICAM,"        << time.nicam() << endl;
    os1 << "Invoke LETKF," << time.letkf_launch() << endl;
    os1 << "LETKF,"        << time.letkf() << endl;
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
    size_t local_count = dps.size();
    size_t total_count;
    MPI_Allreduce(&local_count, &total_count, 1, MPI_LONG, MPI_SUM,
		  env.mpi_comm);
    assert(total_count == (size_t)(kNumEnsemble * kNumCell));
#endif
    time_.letkf_start = gettime(env);

    // alltoall the input data
    int letkf_nprocs;
    int letkf_rank;
    MPI_Comm_size(env.mpi_comm, &letkf_nprocs);
    MPI_Comm_rank(env.mpi_comm, &letkf_rank);
    assert(dps.size() == kNumCell);
    // setup send counts
    int *send_cnts = new int[letkf_nprocs];
    int *sdispls   = new int[letkf_nprocs];
    int each_count = (int)kNumCell / letkf_nprocs;
    int each_rem   = (int)kNumCell % letkf_nprocs;
    for (int i = 0; i < letkf_nprocs; i++) {
      send_cnts[i] =
	(each_count + ((i < each_rem)? 1 : 0)) * (int)kElementCount;
      if (i == 0) {
	sdispls[i] = 0;
      } else {
	sdispls[i] = sdispls[i-1] + send_cnts[i-1];
      }
    }
    assert(sdispls[letkf_nprocs-1] + send_cnts[letkf_nprocs-1] ==
	   (int)(kNumCell * kElementCount));
    // setup send buf
    int sndcnt  = (int)(dps.size() * kElementCount);
    int *sndbuf = new int[sndcnt];
    size_t sndbuf_idx = 0;
    for (vector<DataPack>::iterator itr = dps.begin(); itr != dps.end();
	 itr++) {
      int *data = (int*)itr->data()->value();
      size_t dsize = itr->data()->size();
      memcpy(sndbuf + sndbuf_idx, data, dsize);
      sndbuf_idx += (dsize / sizeof(int));
    }
    // setup receive counts and buf
    int *recv_cnts = new int[letkf_nprocs];
    int *rdispls   = new int[letkf_nprocs];
    int rcvbuf_siz = 0;
    for (int i = 0; i < letkf_nprocs; i++) {
      recv_cnts[i] =
	(each_count + ((letkf_rank < each_rem)? 1 : 0)) * (int)kElementCount;
      rcvbuf_siz += recv_cnts[i];
      if (i == 0) {
	rdispls[i] = 0;
      } else {
	rdispls[i] = rdispls[i-1] + recv_cnts[i-1];
      }
    }
    assert(rcvbuf_siz ==
	   (each_count + ((letkf_rank < each_rem)? 1 : 0)) *
	   (int)(kNumEnsemble * kElementCount));
    int *rcvbuf = new int[rcvbuf_siz];
    MPI_Alltoallv(sndbuf, send_cnts, sdispls, MPI_INT,
		  rcvbuf, recv_cnts, rdispls, MPI_INT, env.mpi_comm);
    delete[] sndbuf;

    // computation: just decrease the value
    for (size_t i = 0; i < (size_t)rcvbuf_siz; i++) {
      rcvbuf[i] -= 1;
    }

    // alltoall the output data
    sndbuf = new int[sndcnt];
    MPI_Alltoallv(rcvbuf, recv_cnts, rdispls, MPI_INT,
		  sndbuf, send_cnts, sdispls, MPI_INT, env.mpi_comm);
    delete[] send_cnts;
    delete[] sdispls;
    delete[] recv_cnts;
    delete[] rdispls;

    // emit data
    int *p = sndbuf;
    for (vector<DataPack>::iterator itr = dps.begin(); itr != dps.end();
	 itr++) {
      Data data(p, itr->data()->size());
      p += (itr->data()->size() / sizeof(int));
      outds->add(itr->key(), data);
    }
    delete[] sndbuf;

    time_.letkf_finish = gettime(env);
    return 0;
  }
};

void run_letkf(DataStore* inds, DataStore* outds, Time& time)
{
  PseudoLETKF mapper(time);
  time.letkf_invoke = gettime();
  View view(kDimEnsembleData);
  bool view_flag[3] = {false, true, false};
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
