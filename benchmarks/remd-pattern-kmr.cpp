/// \file
/// Benchmark program for REMD data access pattern for MPI
///
/// It requires KMR backend, but does not use kmrnext.
/// It assumes that the number of static MPI processes running this program
/// is kNumEnsemble and dynamic MPI processes should be bigger than
/// kNumEnsemble x kNumProc.  It spawns at most once on each process.
#include <iostream>
#include <sstream>
#include <cassert>
#include <ctime>
#include <cstdlib>
#include <unistd.h>
#include <mpi.h>
#include <kmr.h>
#include "../config.hpp"

using namespace std;

const size_t kNumIteration = 5;

#if DEBUG
const size_t kNumEnsemble    = 2;
const size_t kNumProc        = 8;
const size_t kNumData        = 6;
const size_t kElementCount   = 1;    // a double value
const unsigned int kTimeMD   = 1000; // msec
const unsigned int kTimeEX   = 1000; // msec
#else
const size_t kNumEnsemble    = 128;
const size_t kNumProc        = 8;
const size_t kNumData        = 6;
const size_t kElementCount   = 1;
const unsigned int kTimeMD   = 1000; // msec
const unsigned int kTimeEX   = 1000; // msec
#endif

const bool kPrint = true;

int rank             = 0;
int nprocs           = 1;
int universe_size    = 1;
size_t md_task_count = 0;
size_t md_data_count = 0;
int spawn_gap_msec   = 0;

struct Time {
  double loop_start;
  double loop_finish;

  double md_start;
  double md_finish;
  double shuffle_start;
  double shuffle_finish;
  double ex_start;
  double ex_finish;

  double loop() {
#ifdef SYSENV_K_FX
    return (loop_finish - loop_start) / 10E9 - 2;  // -2 for sleep
#else
    return (loop_finish - loop_start) / 10E9;
#endif
  }
  double md() {
    return (md_finish - md_start) / 10E9;
  }
  double shuffle() {
    return (shuffle_finish - shuffle_start) / 10E9;
  }
  double ex() {
    return (ex_finish - ex_start) / 10E9;
  }
};

void setup();
void workflow(string prog_name);
void spawn(string prog_name, string command, int procs, MPI_Comm *icomm);
void disconnect(MPI_Comm *icomm);
void spawn_md(string prog_name, Time& time);
void shuffle(Time& time);
void run_ex(Time& time);
void task_md();
double* allocate_data(size_t count);
double* load_md_data();
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
  string prog_name(argv[0]);
  string subcommand(argv[argc - 1]);
  if (subcommand == "md") {
    task_md();
  } else {
    workflow(prog_name);
  }
  MPI_Finalize();
  return 0;
}

void setup() {
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  int *usizep;
  int uflag;
  MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_UNIVERSE_SIZE, &usizep, &uflag);
  universe_size = *usizep;

  int div_md = (int)kNumEnsemble / nprocs;
  int rem_md = (int)kNumEnsemble % nprocs;
  md_task_count = (rank < rem_md)? div_md + 1 : div_md;
  md_data_count = kNumData * kElementCount;

#if 0
  long spawn_gap_param[2] = {500, 1000};
  // long spawn_gap_param[2] = {1000, 10000};  // KMR original
  int usz = 0;
  unsigned int v = (unsigned int)universe_size;
  while (v > 0) {
    v = (v >> 1);
    usz++;
  }
  spawn_gap_msec = (int)(((spawn_gap_param[1] * usz) / 10)
  			 + spawn_gap_param[0]);
#else
  spawn_gap_msec = 1000;
#endif
}

void workflow(string prog_name) {
#if DEBUG
  assert(nprocs == kNumEnsemble);
  assert(universe_size >= (int)(kNumEnsemble * (1 + kNumProc)));
#endif
  if (rank == 0) {
    cout << "PARAMETERS" << endl;
    cout << "  Universe size     : " << universe_size << endl;
    cout << "  Static proc size  : " << nprocs << endl;
#ifdef SYSENV_K_FX
    cout << "  Spawn gap in msec : " << spawn_gap_msec << endl;
#endif
    cout << endl;
  }

  for (size_t i = 0; i < kNumIteration; i++) {
    Time time;
    time.loop_start = gettime(MPI_COMM_WORLD);

    // run pseudo-MD
    spawn_md(prog_name, time);
#ifdef SYSENV_K_FX
    sleep(2);  // for process release
#endif

    // shuffle
    shuffle(time);

    // run pseudo-EX
    run_ex(time);

    time.loop_finish = gettime(MPI_COMM_WORLD);
    ostringstream os1;
    os1 << "Iteration[" << i << "] ends in " << time.loop() << " sec."
	<< endl;
    os1 << "  MD takes " << time.md() << " sec." << endl;
    os1 << "  Shuffle takes " << time.shuffle() << " sec." << endl;
    os1 << "  EX takes " << time.ex() << " sec." << endl;
    print_line(os1);
#ifdef SYSENV_K_FX
    sleep(2);  // for process release
#endif
  }
}

void spawn(string prog_name, string command, int procs, MPI_Comm *icomm) {
  char *cmd = (char*)prog_name.c_str();
  char *argv[2];
  argv[0] = (char*)command.c_str();
  argv[1] = NULL;
  int *aoe = (int*)calloc(procs, sizeof(int));
  int cc = MPI_Comm_spawn(cmd, argv, procs, MPI_INFO_NULL, 0, MPI_COMM_SELF,
			  icomm, aoe);
  assert(cc == MPI_SUCCESS);
  free(aoe);
}

void disconnect(MPI_Comm *icomm) {
  MPI_Comm_disconnect(icomm);
#ifdef SYSENV_K_FX
  // wait for process release on K
  usleep((useconds_t)(spawn_gap_msec * 1000));
#endif
}

void spawn_md(string prog_name, Time& time) {
  assert(md_task_count == 1);

  time.md_start = gettime(MPI_COMM_WORLD);
  string prog_md("md");
  MPI_Comm icomm;
  spawn(prog_name, prog_md, kNumProc, &icomm);

  int cc;
  double *sbuf = load_md_data();
  cc = MPI_Send(sbuf, (int)md_data_count, MPI_DOUBLE, 0, 1000, icomm);
  assert(cc == MPI_SUCCESS);
  free(sbuf);

  double *rbuf = allocate_data(md_data_count);
  MPI_Status stat;
  cc = MPI_Recv(rbuf, (int)md_data_count, MPI_DOUBLE, 0, 1001, icomm, &stat);
  assert(cc == MPI_SUCCESS);
  free(rbuf);

  disconnect(&icomm);
  time.md_finish = gettime(MPI_COMM_WORLD);
}

void shuffle(Time& time) {
  time.shuffle_start = gettime(MPI_COMM_WORLD);
  KMR *mr = kmr_create_context(MPI_COMM_WORLD, MPI_INFO_NULL, NULL);
  KMR_KVS *kvs0 = kmr_create_kvs(mr, KMR_KV_INTEGER, KMR_KV_OPAQUE);
  double *sbuf = load_md_data();
  struct kmr_kv_box kv;
  kv.klen = (int)sizeof(long);
  kv.vlen = (int)(md_data_count * sizeof(double));
  kv.k.i = 0;
  kv.v.p = (char*)sbuf;
  int cc = kmr_add_kv(kvs0, kv);
  assert(cc == MPI_SUCCESS);
  kmr_add_kv_done(kvs0);

  KMR_KVS *kvs1 = kmr_create_kvs(mr, KMR_KV_INTEGER, KMR_KV_OPAQUE);
  cc = kmr_shuffle(kvs0, kvs1, kmr_noopt);
  assert(cc == MPI_SUCCESS);
  KMR_KVS *kvs2 = kmr_create_kvs(mr, KMR_KV_INTEGER, KMR_KV_OPAQUE);
  cc = kmr_replicate(kvs1, kvs2, kmr_noopt);
  kmr_free_kvs(kvs2);
  free(sbuf);
  time.shuffle_finish = gettime(MPI_COMM_WORLD);
}

void run_ex(Time& time) {
  time.ex_start = gettime(MPI_COMM_WORLD);
  usleep(kTimeEX * 1000);
  time.ex_finish = gettime(MPI_COMM_WORLD);
}

void task_md() {
  assert(nprocs == kNumProc);
  MPI_Comm pcomm;
  MPI_Comm_get_parent(&pcomm);

  if (rank == 0) {
    // receive data from the spawner on rank 0
    double *buf = allocate_data(md_data_count);
    MPI_Status stat;
    int cc = MPI_Recv(buf, (int)md_data_count, MPI_DOUBLE, 0, 1000, pcomm,
		      &stat);
    assert(cc == MPI_SUCCESS);
    free(buf);
  }

  usleep(kTimeMD * 1000);

  if (rank == 0) {
    double *buf = load_md_data();
    // send data to the spawner from rank 0
    int cc = MPI_Send(buf, (int)md_data_count, MPI_DOUBLE, 0, 1001, pcomm);
    assert(cc == MPI_SUCCESS);
    free(buf);
  }

  MPI_Comm_disconnect(&pcomm);
}

double* allocate_data(size_t count) {
  return (double*)calloc(count, sizeof(double));
}

double* load_md_data() {
  double *buf = allocate_data(md_data_count);
  for (size_t i = 0; i < kNumData; i++) {
    size_t idxi = i * kElementCount;
    for (size_t j = 0; j < kElementCount; j++) {
      size_t idx = idxi + j;
      buf[idx] = (int)i;
    }
  }
  return buf;
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
