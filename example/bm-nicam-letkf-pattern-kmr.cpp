/// \file
/// Benchmark program for NICAM-LETKF data access pattern for MPI
///
/// It requires KMR, but does not use kmrnext.
/// It assumes that the number of static MPI processes running this program
/// is kNumEnsemble and dynamic MPI processes should be bigger than
/// kNumEnsemble x kNumRegion.  It spawns at most once on each process at
/// each task, that is, it spawns # of processes times at each task in total.
#include <iostream>
#include <sstream>
#include <cassert>
#include <ctime>
#include <cstdlib>
#include <unistd.h>
#include <mpi.h>
#include <kmr.h>
#include "../config.hpp"

/// If SPAWN_ONCE is 1, each static process spawns only once for LETKF.
/// In that case, the number of spawned MPI processes is same as that of
/// NICAM and the LETKF task performs scatter and gather data at the head
/// and tha end of execution.
///
/// If SPAWN_ONCE is 0, LETKF program is spawned for kNumRegion x kNumCell
/// times and each program runs with 1 MPI process.
#define SPAWN_ONCE 1

using namespace std;

const size_t kNumIteration = 5;

#if DEBUG
const size_t kNumEnsemble	= 2;
const size_t kNumRegion		= 10;
const size_t kNumCell		= 10;
const size_t kElementCount	= 2;
const unsigned int kTimeNICAM	= 1000; // msec
const unsigned int kTimeLETKF	= 1000; // msec
#else
const size_t kNumEnsemble       = 64;
const size_t kNumRegion         = 10;
const size_t kNumCell           = 1156;
// Assume that each lattice has 6160 KB of data (6160 = 1540 * 4)
const size_t kElementCount      = 1540;
const unsigned int kTimeNICAM   = 50000; // msec
const unsigned int kTimeLETKF   = 50;    // msec
#endif

const bool kPrint = true;

int rank                = 0;
int nprocs              = 1;
int universe_size       = 1;
size_t nicam_task_count = 0;
size_t letkf_task_count = 0;
size_t nicam_data_count = 0;
size_t letkf_data_count = 0;
int spawn_gap_msec      = 0;

struct Time {
  double loop_start;
  double loop_finish;

  double nicam_start;
  double nicam_finish;
  double shuffle_start;
  double shuffle_finish;
  double letkf_start;
  double letkf_finish;

  double loop() {
    return (loop_finish - loop_start) / 10E9;
  }
  double nicam() {
    return (nicam_finish - nicam_start) / 10E9;
  }
  double shuffle() {
    return (shuffle_finish - shuffle_start) / 10E9;
  }
  double letkf() {
    return (letkf_finish - letkf_start) / 10E9;
  }
};

void setup();
void workflow(string prog_name);
void spawn(string prog_name, string command, int procs, MPI_Comm *icomm);
void disconnect(MPI_Comm *icomm);
void spawn_nicam(string prog_name, Time& time);
void shuffle(Time& time);
void spawn_letkf(string prog_name, Time& time);
void task_nicam();
void task_letkf();
int* allocate_data(size_t count);
int* load_nicam_data();
int* load_letkf_data();
void print_data(int* data);
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
  if (subcommand == "nicam") {
    task_nicam();
  } else if (subcommand == "letkf") {
    task_letkf();
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

  int div_nicam = (int)kNumEnsemble / nprocs;
  int rem_nicam = (int)kNumEnsemble % nprocs;
  nicam_task_count = (rank < rem_nicam)? div_nicam + 1 : div_nicam;
  int div_letkf = (int)((kNumRegion * kNumCell) / nprocs);
  int rem_letkf = (int)((kNumRegion * kNumCell) % nprocs);
  letkf_task_count = (rank < rem_letkf)? div_letkf + 1 : div_letkf;

  nicam_data_count = kNumRegion * kNumCell * kElementCount;
  letkf_data_count = kNumEnsemble * kElementCount;

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
  assert(universe_size >= (int)(kNumEnsemble * (1 + kNumRegion)));
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

    // run pseudo-NICAM
    spawn_nicam(prog_name, time);
#ifdef SYSENV_K_FX
    sleep(2);  // for process release
#endif

    // shuffle
    shuffle(time);

    // run pseudo-LETKF
    spawn_letkf(prog_name, time);

    time.loop_finish = gettime(MPI_COMM_WORLD);
    ostringstream os1;
    os1 << "Iteration[" << i << "] ends in " << (time.loop() - 2) << " sec."
	<< endl;
    os1 << "  NICAM takes " << time.nicam() << " sec." << endl;
    os1 << "  Shuffle takes " << time.shuffle() << " sec." << endl;
    os1 << "  LETKF takes " << time.letkf() << " sec." << endl;
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

void spawn_nicam(string prog_name, Time& time) {
  assert(nicam_task_count == 1);

  time.nicam_start = gettime(MPI_COMM_WORLD);
  string prog_nicam("nicam");
  MPI_Comm icomm;
  spawn(prog_name, prog_nicam, kNumRegion, &icomm);

  int cc;
  int *sbuf = load_nicam_data();
  cc = MPI_Send(sbuf, (int)nicam_data_count, MPI_INT, 0, 1000, icomm);
  assert(cc == MPI_SUCCESS);
  free(sbuf);

  int *rbuf = allocate_data(nicam_data_count);
  MPI_Status stat;
  cc = MPI_Recv(rbuf, (int)nicam_data_count, MPI_INT, 0, 1001, icomm, &stat);
  assert(cc == MPI_SUCCESS);
  free(rbuf);

  disconnect(&icomm);
  time.nicam_finish = gettime(MPI_COMM_WORLD);
}

void shuffle(Time& time) {
  time.shuffle_start = gettime(MPI_COMM_WORLD);
  KMR *mr = kmr_create_context(MPI_COMM_WORLD, MPI_INFO_NULL, NULL);
  KMR_KVS *kvs0 = kmr_create_kvs(mr, KMR_KV_INTEGER, KMR_KV_OPAQUE);
  int *sbuf = load_nicam_data();
  int cc;
  for (size_t i = 0; i < kNumRegion; i++) {
    for (size_t j = 0; j < kNumCell; j++) {
      size_t idx = i * kNumCell * kElementCount + j * kElementCount;
      struct kmr_kv_box kv;
      kv.klen = (int)sizeof(long);
      kv.vlen = (int)(kElementCount * sizeof(int));
      kv.k.i = idx;
      kv.v.p = (char*)&sbuf[idx];
      cc = kmr_add_kv(kvs0, kv);
      assert(cc == MPI_SUCCESS);
    }
  }
  kmr_add_kv_done(kvs0);

  KMR_KVS *kvs1 = kmr_create_kvs(mr, KMR_KV_INTEGER, KMR_KV_OPAQUE);
  cc = kmr_shuffle(kvs0, kvs1, kmr_noopt);
  assert(cc == MPI_SUCCESS);
  kmr_free_kvs(kvs1);
  free(sbuf);
  time.shuffle_finish = gettime(MPI_COMM_WORLD);
}

struct LETKFProc {
  MPI_Comm icomm;
  int *rbuf;
  MPI_Request rreq;
  bool done;
};

static void wait_letkf_done(LETKFProc *letkf) {
  MPI_Status stat;
  MPI_Wait(&letkf->rreq, &stat);
  free(letkf->rbuf);
  disconnect(&letkf->icomm);
  letkf->done = true;
}

void spawn_letkf(string prog_name, Time& time) {
  time.letkf_start = gettime(MPI_COMM_WORLD);
  string prog_letkf("letkf");
  size_t spawn_count, msg_size;
  int running = 0, nmaxspawn, proc_count;
#if SPAWN_ONCE
  spawn_count = 1;
  nmaxspawn  = 1;
  proc_count = kNumRegion;
  msg_size = letkf_data_count * letkf_task_count;
#else
  spawn_count = letkf_task_count;
  nmaxspawn  = (universe_size - nprocs) / nprocs;
  proc_count = 1;
  msg_size = letkf_data_count;
#endif
  LETKFProc *letkfs = (LETKFProc*)calloc(spawn_count, sizeof(LETKFProc));

  // Spawn all
  for (size_t i = 0; i < spawn_count; i++) {
    if (running == nmaxspawn) {
      // Wait until a dynamic process is freed.
      for (size_t j = 0; j < i; j++) {
	if (letkfs[j].done) { continue; }
	wait_letkf_done(&letkfs[j]);
	running -= 1;
      }
      assert(running < nmaxspawn);
    }

    letkfs[i].done = false;
    spawn(prog_name, prog_letkf, proc_count, &(letkfs[i].icomm));
    // Send data to a spawned process.
    int cc;
#if SPAWN_ONCE
    cc = MPI_Send(&msg_size, 1, MPI_UNSIGNED_LONG, 0, 1000, letkfs[i].icomm);
    assert(cc == MPI_SUCCESS);
#endif
    int *sbuf = load_letkf_data();
    cc = MPI_Send(sbuf, (int)msg_size, MPI_INT, 0, 1000,  letkfs[i].icomm);
    assert(cc == MPI_SUCCESS);
    free(sbuf);

    // Receive data from a spawned process.
    letkfs[i].rbuf = allocate_data(msg_size);
    cc = MPI_Irecv(letkfs[i].rbuf, (int)msg_size, MPI_INT, 0, 1001,
     		   letkfs[i].icomm, &(letkfs[i].rreq));
    assert(cc == MPI_SUCCESS);

    running += 1;
  }

  // Wait all
  for (size_t i = 0; i < spawn_count; i++) {
    if (letkfs[i].done) { continue; }
    wait_letkf_done(&letkfs[i]);
  }

  free(letkfs);
  time.letkf_finish = gettime(MPI_COMM_WORLD);
}

void task_nicam() {
  assert(nprocs == kNumRegion);
  MPI_Comm pcomm;
  MPI_Comm_get_parent(&pcomm);

  if (rank == 0) {
    // receive data from the spawner on rank 0
    int *buf = allocate_data(nicam_data_count);
    MPI_Status stat;
    int cc = MPI_Recv(buf, (int)nicam_data_count, MPI_INT, 0, 1000, pcomm,
		      &stat);
    assert(cc == MPI_SUCCESS);
    free(buf);
  }

  usleep(kTimeNICAM * 1000);

  if (rank == 0) {
    int *buf = load_nicam_data();
    // send data to the spawner from rank 0
    int cc = MPI_Send(buf, (int)nicam_data_count, MPI_INT, 0, 1001, pcomm);
    assert(cc == MPI_SUCCESS);
    free(buf);
  }

  MPI_Comm_disconnect(&pcomm);
}

void task_letkf() {
  MPI_Comm pcomm;
  MPI_Comm_get_parent(&pcomm);
  size_t msg_size;
#if SPAWN_ONCE
  assert(nprocs == kNumRegion);
#else
  assert(nprocs == 1);
  msg_size = letkf_data_count;
#endif
  int cc;
  int *buf = NULL;

  if (rank == 0) {
    // receive data from the spawner on rank 0
    MPI_Status stat;
#if SPAWN_ONCE
    cc = MPI_Recv(&msg_size, 1, MPI_UNSIGNED_LONG, 0, 1000, pcomm, &stat);
    assert(cc == MPI_SUCCESS);
#endif
    buf = allocate_data(msg_size);
    cc = MPI_Recv(buf, (int)msg_size, MPI_INT, 0, 1000, pcomm, &stat);
    assert(cc == MPI_SUCCESS);
  }

#if SPAWN_ONCE
  MPI_Bcast(&msg_size, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
  int each_size = (int)msg_size / nprocs;
  int *ebuf = (int*)calloc(each_size, sizeof(int));
  cc = MPI_Scatter(buf, each_size, MPI_INT, ebuf, each_size, MPI_INT, 0,
		   MPI_COMM_WORLD);
#endif

  usleep(kTimeLETKF * 1000);

#if SPAWN_ONCE
  int *rbuf;
  if (rank == 0) {
    rbuf = (int*)calloc(each_size * nprocs, sizeof(int));
  }
  cc = MPI_Gather(ebuf, each_size, MPI_INT, rbuf, each_size, MPI_INT, 0,
		  MPI_COMM_WORLD);
  free(ebuf);
  if (rank == 0) {
    free(rbuf);
  }
#endif

  if (rank == 0) {
    // send data to the spawner from rank 0
    cc = MPI_Send(buf, (int)msg_size, MPI_INT, 0, 1001, pcomm);
    assert(cc == MPI_SUCCESS);
    free(buf);
  }

  MPI_Comm_disconnect(&pcomm);
}

int* allocate_data(size_t count) {
  return (int*)calloc(count, sizeof(int));
}

int* load_nicam_data() {
  int *buf = allocate_data(nicam_data_count);
  for (size_t i = 0; i < kNumRegion; i++) {
    size_t idxi = i * kNumCell * kElementCount;
    for (size_t j = 0; j < kNumCell; j++) {
      size_t idxj = j * kElementCount;
      for (size_t k = 0; k < kElementCount; k++) {
	size_t idx = idxi + idxj + k;
	buf[idx] = (int)i;
      }
    }
  }
  return buf;
}

int* load_letkf_data() {
#if SPAWN_ONCE
  int *buf = allocate_data(letkf_data_count * letkf_task_count);
  for (size_t h = 0; h < letkf_task_count; h++) {
    size_t idxh = h * letkf_data_count;
    for (size_t i = 0; i < kNumEnsemble; i++) {
      size_t idxi = i * kElementCount;
      for (size_t k = 0; k < kElementCount; k++) {
	size_t idx = idxh + idxi + k;
	buf[idx] = (int)i;
      }
    }
  }
#else
  int *buf = allocate_data(letkf_data_count);
  for (size_t i = 0; i < kNumEnsemble; i++) {
    size_t idxi = i * kElementCount;
    for (size_t k = 0; k < kElementCount; k++) {
      size_t idx = idxi + k;
      buf[idx] = (int)i;
    }
  }
#endif
  return buf;
}

void print_data(int* data) {
  const int p2ptag = 1000;
  for (int i = 0; i < nprocs; i++) {
    if (i == 0) {
      // just print on rank 0
      if (rank != 0) { continue; }
      for (size_t j = 0; j < kNumRegion; j++) {
	for (size_t k = 0; k < kNumCell; k++) {
	  for (size_t l = 0; l < kElementCount; l++) {
	    cout << data[j * kNumCell * kElementCount + k * kElementCount + l]
		 << " ";
	  }
	  cout << endl;
	}
      }
      cout << endl;
    } else {
      // send data to rank 0 and print on rank 0
      if (rank == 0) {
	int *buf = (int*)calloc(kElementCount, sizeof(int));
	MPI_Status stat;
	MPI_Recv(buf, kElementCount, MPI_INT, i, p2ptag, MPI_COMM_WORLD,
		 &stat);
	for (size_t j = 0; j < kNumRegion; j++) {
	  for (size_t k = 0; k < kNumCell; k++) {
	    for (size_t l = 0; l < kElementCount; l++) {
	      cout << data[j * kNumCell * kElementCount + k * kElementCount + l]
		   << " ";
	    }
	    cout << endl;
	  }
	}
	cout << endl;
      } else if (rank == i) {
	MPI_Send(data, kElementCount, MPI_INT, 0, p2ptag, MPI_COMM_WORLD);
      }
    }
  }
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
