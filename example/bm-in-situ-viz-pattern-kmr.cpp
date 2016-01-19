/// \file
/// Benchmark program for visualization data access pattern for MPI
///
/// It only requires MPI.
/// The number of process can divide kX x kY x kZ.
#include <iostream>
#include <sstream>
#include <cassert>
#include <ctime>
#include <cstdlib>
#include <unistd.h>
#include <mpi.h>

using namespace std;

const size_t kDimSpace = 3;

#if DEBUG
const size_t kX             = 2;
const size_t kY             = 2;
const size_t kZ	            = 2;
const size_t kDataCount	    = 2;
const unsigned int kTimeSim = 1; // sec
const unsigned int kTimeViz = 1; // sec
#else
const size_t kX             = 1000;
const size_t kY             = 1000;
const size_t kZ             = 1000;
// Assume that each lattice has 1000 KB of data (1000 = 250 * 4)
const size_t kDataCount     = 250;
const unsigned int kTimeSim = 60;  // sec
const unsigned int kTimeViz = 1;   // sec
#endif

const size_t kSpaceSizes[kDimSpace] = {kX, kY, kZ};

const bool kPrint = true;

int rank                = 0;
int nprocs              = 1;
size_t total_data_count = 0;

struct Time {
  double main_start;
  double main_finish;

  double sim_start;
  double sim_finish;
  double viz_start;
  double viz_finish;

  double whole() {
    return (main_finish - main_start) / 10E9;
  }
  double sim() {
    return (sim_finish - sim_start) / 10E9;
  }
  double viz() {
    return (viz_finish - viz_start) / 10E9;
  }
};

void setup();
int* load_data();
int* allocate_data();
void run_simulation(int* in, int* out, Time& time);
void run_viz(int* in, int* out, Time& time);
void print_line(string& str);
void print_line(ostringstream& os);
double gettime(MPI_Comm comm);

//////////////////////////////////////////////////////////////////////////////
// Main starts from here.
//////////////////////////////////////////////////////////////////////////////
int
main(int argc, char **argv)
{
  int thlv;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &thlv);
  setup();

  Time time;
  time.main_start = gettime(MPI_COMM_WORLD);

  // load data
  int *data0 = load_data();

  // run pseudo-Simulation
  int *data1 = allocate_data();
  run_simulation(data0, data1, time);
  free(data0);

  // run pseudo-Visualization
  run_viz(data1, NULL, time);
  free(data1);

  time.main_finish = gettime(MPI_COMM_WORLD);
  ostringstream os;
  os << "Total time: " << time.whole() << " sec." << endl;
  os << "  Simulation consumes " << time.sim() << " sec." << endl;
  os << "  Visualization consumes " << time.viz() << " sec." << endl;
  print_line(os);

  MPI_Finalize();
  return 0;
}


void setup() {
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
#if DEBUG
  assert((kX * kY * kZ) % nprocs == 0);
#endif
  total_data_count = kDataCount;
  for (size_t i = 0; i < kDimSpace; i++) {
    total_data_count *= kSpaceSizes[i];
  }
}

int* allocate_data() {
  int* data_ = NULL;
  if (rank == 0) {
    data_ = (int*)calloc(total_data_count, sizeof(int));
  }
  return data_;
}

int* load_data() {
  int* data_ = NULL;
  if (rank == 0) {
    data_ = (int*)calloc(total_data_count, sizeof(int));
    for (size_t i = 0; i < total_data_count; i++) {
      data_[i] = rank + 1;
    }
  }
  return data_;
}

void run_simulation(int* in, int* out, Time& time)
{
  time.sim_start = gettime(MPI_COMM_WORLD);

  int send_cnt = (int)total_data_count / nprocs;
#if DEBUG
  assert(send_cnt == (int)(kDataCount * ((kX * kY * kZ) / (size_t)nprocs)));
#endif
  int *rbuf = (int*)calloc(send_cnt, sizeof(int));
  MPI_Scatter(in, send_cnt, MPI_INT, rbuf, send_cnt, MPI_INT,
	      0, MPI_COMM_WORLD);
#if DEBUG
  for (size_t i = 0; i < (size_t)send_cnt; i++) {
    rbuf[i] += 1;
  }
#endif
  sleep(kTimeSim);
  MPI_Gather(rbuf, send_cnt, MPI_INT, out, send_cnt, MPI_INT,
	     0, MPI_COMM_WORLD);

  time.sim_finish = gettime(MPI_COMM_WORLD);
}

void run_viz(int* in, int* out, Time& time)
{
  time.viz_start = gettime(MPI_COMM_WORLD);

  int send_cnt = (int)total_data_count / nprocs;
#if DEBUG
  assert(send_cnt == (int)(kDataCount * ((kX * kY * kZ) / (size_t)nprocs)));
#endif
  int *rbuf = (int*)calloc(send_cnt, sizeof(int));
  MPI_Scatter(in, send_cnt, MPI_INT, rbuf, send_cnt, MPI_INT,
	      0, MPI_COMM_WORLD);
#if DEBUG
  for (size_t i = 0; i < (size_t)send_cnt; i++) {
    rbuf[i] += 1;
  }
#endif
  sleep(kTimeViz);

  time.viz_finish = gettime(MPI_COMM_WORLD);
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

double gettime(MPI_Comm comm) {
  MPI_Barrier(comm);
  struct timespec ts;
  clock_gettime(CLOCK_REALTIME, &ts);
  return ((double) ts.tv_sec) * 10E9 + ((double) ts.tv_nsec);
}
