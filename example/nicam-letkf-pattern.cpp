#include <iostream>
#include <sstream>
#include <cassert>
#include <ctime>
#include "kmrnext.hpp"

using namespace std;
using namespace kmrnext;

int rank = 0;
int nprocs = 1;

// If true, LETKF works on regions. If false, it works on each cell.
//
// When LETKF works on regions, LETKF uses multiple nodes.
// When LETKF works on cells, LETKF uses single node.
const bool kLETKFRegion = false;

const size_t kNumIteration = 10;

const size_t kDimEnsembleData = 3;

#if DEBUG
const size_t kNumEnsemble  = 2;
const size_t kNumRegion    = 10;
const size_t kNumCell      = 10;
const size_t kElementCount = 2;
#else
const size_t kNumEnsemble  = 64;
const size_t kNumRegion    = 10;
const size_t kNumCell      = 1156;
// Assume that each grid point has 6160 bytes of data (6160 = 1540 * 4)
const size_t kElementCount = 1540;
#endif

const size_t kEnsembleDataDimSizes[kDimEnsembleData] =
  {kNumEnsemble, kNumRegion, kNumCell};

const bool kPrint = true;

#if DEBUG
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
#endif

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
void run_nicam(DataStore* ds, Time& time);
void run_letkf(DataStore* ds, Time& time);
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
#ifdef BACKEND_KMR
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
#endif

  DataStore* ds0 = next->create_ds(kDimEnsembleData);
  ds0->set(kEnsembleDataDimSizes);
  load_data(ds0);
#if DEBUG
  DataPrinter dp;
  string str0 = ds0->dump(dp);
  print_line(str0);
#endif

  for (size_t i = 0; i < kNumIteration; i++) {
    ostringstream os0;
    os0 << "Iteration[" << i << "] starts.";
    print_line(os0);
    Time time;
    time.loop_start = gettime();

    // run pseudo-NICAM
    run_nicam(ds0, time);
#if DEBUG
    string str1 = ds0->dump(dp);
    print_line(str1);
#endif

    // run pseudo-LETKF
    run_letkf(ds0, time);
#if DEBUG
    string str2 = ds0->dump(dp);
    print_line(str2);
#endif

    time.loop_finish = gettime();
    ostringstream os1;
    os1 << "Iteration[" << i << "] ends in " << time.loop() << " sec." << endl;
    os1 << "  Invoking NICAM takes " << time.nicam_launch() << " sec." << endl;
    os1 << "    NICAM consumes " << time.nicam() << " sec." << endl;
    os1 << "  Invoking LETKF takes " << time.letkf_launch() << " sec." << endl;
    os1 << "    LETKF consumes " << time.letkf() << " sec." << endl;
    print_line(os1);
  }
  delete ds0;

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
    for (size_t i = 0; i < kElementCount; i++) {
      data_val[i] = (int)y + 1;
    }
    Data data((void*)data_val, sizeof(int) * kElementCount);

    Key key(kDimEnsembleData);
    key.set_dim(0, x);
    key.set_dim(1, y);
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
#ifdef BACKEND_SERIAL
    assert(dps.size() == (size_t)(kNumRegion * kNumCell));
#elif defined BACKEND_KMR
    int nprocs_nicam;
    MPI_Comm_size(env.mpi_comm, &nprocs_nicam);
    int nprocs_calc = calculate_task_nprocs(env.view, env.split,
					    nprocs_nicam);
    assert(nprocs_nicam == nprocs_calc);
    size_t local_count = dps.size();
    size_t total_count;
    MPI_Allreduce(&local_count, &total_count, 1, MPI_LONG, MPI_SUM,
		  env.mpi_comm);
    assert(total_count == (size_t)(kNumRegion * kNumCell));
#endif
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

void run_nicam(DataStore* ds, Time& time)
{
  PseudoNICAM mapper(time);
#ifdef BACKEND_KMR
  View split(kDimEnsembleData);
  bool split_flag[3] = {true, true, false};
  split.set(split_flag);
  ds->set_split(split);
#endif
  time.nicam_invoke = gettime();
  View view(kDimEnsembleData);
  bool view_flag[3] = {true, false, false};
  view.set(view_flag);
  ds->map(mapper, view);
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
#ifdef BACKEND_SERIAL
    if (kLETKFRegion) {
      assert(dps.size() == (size_t)(kNumEnsemble * kNumCell));
    } else {
      assert(dps.size() == (size_t)kNumEnsemble);
    }
#elif defined BACKEND_KMR
    int nprocs_letkf;
    MPI_Comm_size(env.mpi_comm, &nprocs_letkf);
    int nprocs_calc = calculate_task_nprocs(env.view, env.split,
					    nprocs_letkf);
    assert(nprocs_letkf == nprocs_calc);
    size_t local_count = dps.size();
    size_t total_count;
    MPI_Allreduce(&local_count, &total_count, 1, MPI_LONG, MPI_SUM,
		  env.mpi_comm);
    if (kLETKFRegion) {
      assert(total_count == (size_t)(kNumEnsemble * kNumCell));
    } else {
      assert(total_count == (size_t)kNumEnsemble);
    }
#endif
#endif

    time_.letkf_start = gettime(env);
    for (vector<DataPack>::iterator itr = dps.begin(); itr != dps.end();
	 itr++) {
      int *data_new = new int[kElementCount];
      int *data_old = (int*)itr->data()->value();
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

void run_letkf(DataStore* ds, Time& time)
{
  PseudoLETKF mapper(time);
  time.letkf_invoke = gettime();
  View view(kDimEnsembleData);
  if (kLETKFRegion) {
#ifdef BACKEND_KMR
    View split(kDimEnsembleData);
    bool split_flag[3] = {true, true, false};
    split.set(split_flag);
    ds->set_split(split);
#endif
    bool view_flag[3] = {false, true, false};
    view.set(view_flag);
  } else {
#ifdef BACKEND_KMR
    View split(kDimEnsembleData);
    bool split_flag[3] = {false, true, true};
    split.set(split_flag);
    ds->set_split(split);
#endif
    bool view_flag[3] = {false, true, true};
    view.set(view_flag);
  }
  ds->map(mapper, view);
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

int calculate_task_nprocs(View& view, View& alc_view, int given_nprocs) {
  int total_nprocs = 1;
  int nprocs_calc = 1;
  for (size_t i = 0; i < view.size(); i++) {
    if (!view.dim(i) && alc_view.dim(i)) {
      nprocs_calc *= (int)kEnsembleDataDimSizes[i];
    }
    if (alc_view.dim(i)) {
      total_nprocs *= (int)kEnsembleDataDimSizes[i];
    }
  }
  if (nprocs < total_nprocs) {
    nprocs_calc = given_nprocs;
  }
  return nprocs_calc;
}
