#include <iostream>
#include <sstream>
#include <cassert>
#include <ctime>
#include "kmrnext.hpp"

using namespace std;
using namespace kmrnext;

int rank = 0;
int nprocs = 1;

const size_t kDimSpace = 3;

#if DEBUG
const size_t kX		= 2;
const size_t kY		= 2;
const size_t kZ		= 2;
const size_t kDataCount = 2;
#else
const size_t kX         = 100;
const size_t kY         = 100;
const size_t kZ         = 100;
// Assume that each lattice has 4000 KB of data (4000 = 1000 * 4)
// In total, 4GB of data (kX x kY x kZ x kDataCount)
const size_t kDataCount = 1000;
#endif

const size_t kSpaceSizes[kDimSpace] = {kX, kY, kZ};

const bool kPrint = true;

#if DEBUG
class DataPrinter : public DataPack::Dumper {
public:
  string operator()(DataPack& dp)
  {
    ostringstream os;
    os << dp.key().to_string() << " : ";
    int *data = static_cast<int*>(dp.data()->value());
    for (size_t i = 0; i < kDataCount; i++) {
      os << data[i] << " ";
    }
    os << endl;
    return os.str();
  }
};
#endif

struct Time {
  double main_start;
  double main_finish;

  double sim_start;
  double sim_finish;
  double viz_start;
  double viz_finish;

  double sim_invoke;
  double sim_cleanup;
  double viz_invoke;
  double viz_cleanup;

  double whole() {
    return (main_finish - main_start) / 10E9;
  }
  double sim() {
    return (sim_finish - sim_start) / 10E9;
  }
  double viz() {
    return (viz_finish - viz_start) / 10E9;
  }
  double sim_launch() {
    return (sim_cleanup - sim_invoke) / 10E9;
  }
  double viz_launch() {
    return (viz_cleanup - viz_invoke) / 10E9;
  }
};

void load_data(DataStore* ds);
void run_simulation(DataStore* ds, Time& time);
void run_viz(DataStore* ds, Time& time);
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

  Time time;
  time.main_start = gettime();

  // load data
  DataStore* ds0 = next->create_ds(kDimSpace);
  ds0->set(kSpaceSizes);
  load_data(ds0);
#if DEBUG
  DataPrinter dp;
  string str0 = ds0->dump(dp);
  print_line(str0);
#endif

  // run pseudo-Simulation
  run_simulation(ds0, time);
#if DEBUG
  string str1 = ds0->dump(dp);
  print_line(str1);
#endif

  // run pseudo-Visualization
  run_viz(ds0, time);
  delete ds0;

  time.main_finish = gettime();
  ostringstream os;
  os << "Total time: " << time.whole() << " sec." << endl;
  os << "  Invoking Simulation takes " << time.sim_launch() << " sec." << endl;
  os << "    Simulation consumes " << time.sim() << " sec." << endl;
  os << "  Invoking Visualization takes " << time.viz_launch() << " sec."
     << endl;
  os << "    Visualization consumes " << time.viz() << " sec." << endl;
  print_line(os);

  KMRNext::finalize();
  return 0;
}


class DataLoader : public DataStore::Loader<long> {
public:
  int operator()(DataStore* ds, const long& num)
  {
    size_t x = num / kY;
    size_t y = num % kY;
    int *data_val = new int[kDataCount];
    for (size_t i = 0; i < kDataCount; i++) {
      data_val[i] = static_cast<int>(y) + 1;
    }
    Data data(static_cast<void*>(data_val), sizeof(int) * kDataCount);

    Key key(kDimSpace);
    key.set_dim(0, x);
    key.set_dim(1, y);
    for (size_t i = 0; i < kZ; i++) {
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
  for (size_t i = 0; i < kX; i++) {
    for (size_t j = 0; j < kY; j++) {
      data_srcs.push_back(i * kY + j);
    }
  }
  DataLoader dl;
  ds->load_integers(data_srcs, dl);
}

class PseudoSimulation : public DataStore::Mapper {
  Time& time_;
public:
  PseudoSimulation(Time& time) : time_(time) {};

  int operator()(DataStore* inds, DataStore* outds,
		 Key& key, vector<DataPack>& dps,
		 DataStore::MapEnvironment& env)
  {
#if DEBUG
#ifdef BACKEND_SERIAL
    assert(dps.size() == static_cast<size_t>(kX * kY * kZ));
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
    assert(total_count == static_cast<size_t>(kX * kY * kZ));
#endif
#endif

    time_.sim_start = gettime(env);
    for (vector<DataPack>::iterator itr = dps.begin(); itr != dps.end();
	 itr++) {
      int *data_new = new int[kDataCount];
      int *data_old = static_cast<int*>(itr->data()->value());
      for (size_t i = 0; i < kDataCount; i++) {
	data_new[i] = data_old[i] + 1;
      }
      Data data(data_new, itr->data()->size());
      outds->add(itr->key(), data);
      delete[] data_new;
    }
    time_.sim_finish = gettime(env);
    return 0;
  }
};

void run_simulation(DataStore* ds, Time& time)
{
  PseudoSimulation mapper(time);
  time.sim_invoke = gettime();
#ifdef BACKEND_KMR
  View split(kDimSpace);
  bool split_flag[3] = {true, true, false};
  split.set(split_flag);
  ds->set_split(split);
#endif
  View view(kDimSpace);
  bool view_flag[3] = {false, false, false};
  view.set(view_flag);
  ds->map(mapper, view);
  time.sim_cleanup = gettime();
}

class PseudoVisualization : public DataStore::Mapper {
  Time& time_;
public:
  PseudoVisualization(Time& time) : time_(time) {};

  int operator()(DataStore* inds, DataStore* outds,
		 Key& key, vector<DataPack>& dps,
		 DataStore::MapEnvironment& env)
  {
#if DEBUG
#ifdef BACKEND_SERIAL
    assert(dps.size() == 1);
#elif defined BACKEND_KMR
    int nprocs_viz;
    MPI_Comm_size(env.mpi_comm, &nprocs_viz);
    assert(nprocs_viz == 1);
    assert(dps.size() == 1);
#endif
#endif

    time_.viz_start = gettime(env);
    for (vector<DataPack>::iterator itr = dps.begin(); itr != dps.end();
	 itr++) {
      int *data_new = new int[kDataCount];
      int *data_old = static_cast<int*>(itr->data()->value());
      for (size_t i = 0; i < kDataCount; i++) {
	data_new[i] = data_old[i] - 1;
      }
      Data data(data_new, itr->data()->size());
      delete[] data_new;
    }
    time_.viz_finish = gettime(env);
    return 0;
  }
};

void run_viz(DataStore* ds, Time& time)
{
  PseudoVisualization mapper(time);
  time.viz_invoke = gettime();
#ifdef BACKEND_KMR
  // As task runs on each point, any Split is OK.
#endif
  View view(kDimSpace);
  bool view_flag[3] = {true, true, true};
  view.set(view_flag);
  ds->map(mapper, view);
  time.viz_cleanup = gettime();
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
  return (static_cast<double>(ts.tv_sec) * 10E9 +
	  static_cast<double>(ts.tv_nsec));
}

double gettime(DataStore::MapEnvironment& env) {
#ifdef BACKEND_KMR
  MPI_Barrier(env.mpi_comm);
#endif
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
      nprocs_calc *= static_cast<int>(kSpaceSizes[i]);
    }
    if (alc_view.dim(i)) {
      total_nprocs *= static_cast<int>(kSpaceSizes[i]);
    }
  }
  if (nprocs < total_nprocs) {
    nprocs_calc = given_nprocs;
  }
  return nprocs_calc;
}
