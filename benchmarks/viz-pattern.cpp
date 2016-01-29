/// \file
/// Benchmark program for in-situ visualization data access pattern
///
/// It requires KMR backend.
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

const size_t kDimSpace = 3;
const size_t kDimData  = 1;

#if DEBUG
const size_t kX             = 2;
const size_t kY             = 2;
const size_t kZ             = 2;
const size_t kDataCount     = 2;
const unsigned int kTimeSim = 1000;  // msec
const unsigned int kTimeViz = 1000;  // msec
#else
const size_t kX             = 128;
const size_t kY             = 128;
const size_t kZ             = 128;
// Assume that each point has 2048 Bytes of data (2048 = 512 * 4)
// In total, 4GB of data (kX x kY x kZ x kDataCount)
const size_t kDataCount     = 512;
const unsigned int kTimeSim = 60000; // msec
const unsigned int kTimeViz = 100;   // msec
#endif

const size_t kSpaceSizes[kDimSpace] = {kX, kY, kZ};

const bool kPrint = true;

class DataPrinter : public DataPack::Dumper {
public:
  string operator()(DataPack& dp)
  {
    ostringstream os;
    os << dp.key().to_string() << " : ";
    int *data = (int*)dp.data()->value();
    for (size_t i = 0; i < kDataCount; i++) {
      os << data[i] << " ";
    }
    os << endl;
    return os.str();
  }
};

struct Time {
  double main_start;
  double main_finish;

  double load_start;
  double load_finish;
  double sim_start;
  double sim_finish;
  double viz_start;
  double viz_finish;

  double whole() {
    return (main_finish - main_start) / 10E9;
  }
  double load() {
    return (load_finish - load_start) / 10E9;
  }
  double sim() {
    return (sim_finish - sim_start) / 10E9;
  }
  double viz() {
    return (viz_finish - viz_start) / 10E9;
  }
};

void load_data(DataStore* ds, Time& time);
void run_simulation(DataStore* inds, DataStore* outds, Time& time);
void run_viz(DataStore* inds, DataStore* outds, Time& time);
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
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  Time time;
  time.main_start = gettime();

  // load data
  DataStore* ds0 = next->create_ds(kDimSpace);
  ds0->set(kSpaceSizes);
  load_data(ds0, time);
  DataPrinter dp;

  // run pseudo-Simulation
  DataStore* ds1 = next->create_ds(kDimSpace);
  ds1->set(kSpaceSizes);
  run_simulation(ds0, ds1, time);
  delete ds0;

  // run pseudo-Visualization
  run_viz(ds1, NULL, time);
  delete ds1;

  time.main_finish = gettime();
  ostringstream os;
  os << "Total time: " << time.whole() << " sec." << endl;
  os << "  Loading data consumes " << time.load() << " sec." << endl;
  os << "  Simulation consumes " << time.sim() << " sec." << endl;
  os << "  Visualization consumes " << time.viz() << " sec." << endl;
  print_line(os);

  KMRNext::finalize();
  return 0;
}


class DataLoader : public DataStore::Loader<int> {
public:
  int operator()(DataStore* ds, const int& num)
  {
    Key key(kDimData);
    int *data_val = new int[kDataCount];
    for (size_t i = 0; i < kDataCount; i++) {
      data_val[i] = num;
    }
    Data data((void*)data_val, sizeof(int) * kDataCount);

    for (size_t i = 0; i < kZ; i++) {
      key.set_dim(0, i);
      ds->add(key, data);
    }
    delete[] data_val;
    return 0;
  }
};

void load_data(DataStore* ds, Time& time)
{
  time.load_start = gettime();
  vector<int> data_srcs;
  for (size_t i = 0; i < kX; i++) {
    for (size_t j = 0; j < kY; j++) {
      data_srcs.push_back((int)(j+1));
    }
  }
  DataLoader dl;
  ds->load_array(data_srcs, dl);
  time.load_finish = gettime();
}

class PseudoSimulation : public DataStore::Mapper {
public:
  int operator()(DataStore* inds, DataStore* outds,
		 Key& key, vector<DataPack>& dps,
		 DataStore::MapEnvironment& env)
  {
#if DEBUG
    int nprocs_sim;
    MPI_Comm_size(env.mpi_comm, &nprocs_sim);
    assert(nprocs_sim == nprocs);
    size_t local_count = dps.size();
    size_t total_count;
    MPI_Allreduce(&local_count, &total_count, 1, MPI_LONG, MPI_SUM,
		  env.mpi_comm);
    assert(total_count == (size_t)(kX * kY * kZ));
#endif

    usleep(kTimeSim * 1000);

    for (vector<DataPack>::iterator itr = dps.begin(); itr != dps.end();
	 itr++) {
      int *data_new = new int[kDataCount];
      int *data_old = (int*)itr->data()->value();
      for (size_t i = 0; i < kDataCount; i++) {
	data_new[i] = data_old[i] + 1;
      }
      Data data(data_new, itr->data()->size());
      outds->add(itr->key(), data);
      delete[] data_new;
    }
    return 0;
  }
};

void run_simulation(DataStore* inds, DataStore* outds, Time& time)
{
  PseudoSimulation mapper;
  time.sim_start = gettime();
  View view(kDimSpace);
  bool view_flag[3] = {false, false, false};
  view.set(view_flag);
  inds->map(outds, mapper, view);
  time.sim_finish = gettime();
}

class PseudoVisualization : public DataStore::Mapper {
  int operator()(DataStore* inds, DataStore* outds,
		 Key& key, vector<DataPack>& dps,
		 DataStore::MapEnvironment& env)
  {
#if DEBUG
    int nprocs_viz;
    MPI_Comm_size(env.mpi_comm, &nprocs_viz);
    assert(nprocs_viz == 1);
    assert(dps.size() == 1);
#endif

    usleep(kTimeViz * 1000);

    return 0;
  }
};

void run_viz(DataStore* inds, DataStore* outds, Time& time)
{
  PseudoVisualization mapper;
  time.viz_start = gettime();
  View view(kDimSpace);
  bool view_flag[3] = {true, true, true};
  view.set(view_flag);
  inds->map(outds, mapper, view);
  time.viz_finish = gettime();
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
