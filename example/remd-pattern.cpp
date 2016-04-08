#include <iostream>
#include <sstream>
#include <cassert>
#include <ctime>
#include "kmrnext.hpp"

using namespace std;
using namespace kmrnext;

int rank = 0;
int nprocs = 1;

const size_t kNumIteration = 10;

const size_t kNumDimensions = 3;
#if DEBUG
const size_t kNumEnsemble   = 2;
const size_t kNumProc       = 8;
const size_t kNumData       = 6;
const size_t kElementCount  = 1;  // a double value
#else
const size_t kNumEnsemble   = 128;
const size_t kNumProc       = 8;
const size_t kNumData       = 6;
const size_t kElementCount  = 1;  // a double value
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
void run_md(DataStore* ds, Time& time);
void run_ex(DataStore* ds, Time& time);
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

  DataStore* ds0 = next->create_ds(kNumDimensions);
  ds0->set(kDataStoreSizes);
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

    // run pseudo-MD
    run_md(ds0, time);
#if DEBUG
    string str1 = ds0->dump(dp);
    print_line(str1);
#endif

    // run pseudo-EX
    run_ex(ds0, time);
#if DEBUG
    string str2 = ds0->dump(dp);
    print_line(str2);
#endif

    time.loop_finish = gettime();
    ostringstream os1;
    os1 << "Iteration[" << i << "] ends in " << time.loop() << " sec." << endl;
    os1 << "  Invoking MD takes " << time.md_launch() << " sec." << endl;
    os1 << "    MD consumes " << time.md() << " sec." << endl;
    os1 << "  Invoking EX takes " << time.ex_launch() << " sec." << endl;
    os1 << "    EX consumes " << time.ex() << " sec." << endl;
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
    size_t x = num / kNumProc;
    size_t y = num % kNumProc;
    double *data_val = new double[kElementCount];
    for (size_t i = 0; i < kElementCount; i++) {
      data_val[i] = (double)(y + 1);
    }
    Data data((void*)data_val, sizeof(double) * kElementCount);

    Key key(kNumDimensions);
    key.set_dim(0, x);
    key.set_dim(1, y);
    for (size_t i = 0; i < kNumData; i++) {
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
    for (size_t j = 0; j < kNumProc; j++) {
      data_srcs.push_back((int)(i * kNumProc + j));
    }
  }
  DataLoader dl;
  ds->load_integers(data_srcs, dl);
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
#ifdef BACKEND_SERIAL
    assert(dps.size() == (size_t)(kNumProc * kNumData));
#elif defined BACKEND_KMR
    int nprocs_sim;
    MPI_Comm_size(env.mpi_comm, &nprocs_sim);
    int nprocs_calc = calculate_task_nprocs(env.view, env.split,
					    nprocs_sim);
    assert(nprocs_sim == nprocs_calc);
    size_t local_count = dps.size();
    size_t total_count;
    MPI_Allreduce(&local_count, &total_count, 1, MPI_LONG, MPI_SUM,
		  env.mpi_comm);
    assert(total_count == (size_t)(kNumProc * kNumData));
#endif
#endif

    time_.md_start = gettime(env);
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
    time_.md_finish = gettime(env);
    return 0;
  }
};

void run_md(DataStore* ds, Time& time)
{
  PseudoMD mapper(time);
  time.md_invoke = gettime();
#ifdef BACKEND_KMR
  View split(kNumDimensions);
  bool split_flag[3] = {true, true, false};
  split.set(split_flag);
  ds->set_split(split);
#endif
  View view(kNumDimensions);
  bool view_flag[3] = {true, false, false};
  view.set(view_flag);
  ds->map(mapper, view);
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
#ifdef BACKEND_SERIAL
      assert(dps.size() == (size_t)(kNumEnsemble * kNumData));
#elif defined BACKEND_KMR
    int nprocs_ex;
    MPI_Comm_size(env.mpi_comm, &nprocs_ex);
    int nprocs_calc = calculate_task_nprocs(env.view, env.split,
					    nprocs_ex);
    assert(nprocs_ex == nprocs_calc);
    size_t local_count = dps.size();
    size_t total_count;
    MPI_Allreduce(&local_count, &total_count, 1, MPI_LONG, MPI_SUM,
		  env.mpi_comm);
      assert(total_count == (size_t)(kNumEnsemble * kNumData));
#endif
#endif

    time_.ex_start = gettime(env);
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
    time_.ex_finish = gettime(env);
    return 0;
  }
};

void run_ex(DataStore* ds, Time& time)
{
  PseudoEX mapper(time);
#ifdef BACKEND_KMR
  View split(kNumDimensions);
  bool split_flag[3] = {false, true, false};
  split.set(split_flag);
  ds->set_split(split);
#endif
  time.ex_invoke = gettime();
  View view(kNumDimensions);
  bool view_flag[3] = {false, true, false};
  view.set(view_flag);
  ds->map(mapper, view);
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

int calculate_task_nprocs(View& view, View& alc_view, int given_nprocs) {
  int total_nprocs = 1;
  int nprocs_calc = 1;
  for (size_t i = 0; i < view.size(); i++) {
    if (!view.dim(i) && alc_view.dim(i)) {
      nprocs_calc *= (int)kDataStoreSizes[i];
    }
    if (alc_view.dim(i)) {
      total_nprocs *= (int)kDataStoreSizes[i];
    }
  }
  if (nprocs < total_nprocs) {
    nprocs_calc = given_nprocs;
  }
  return nprocs_calc;
}
