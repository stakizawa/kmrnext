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

using namespace std;

const size_t kNumIteration = 20;

#if DEBUG
const size_t kNumEnsemble	= 2;
const size_t kNumRegion		= 10;
const size_t kNumCell		= 10;
const size_t kElementCount	= 2;
const unsigned int kTimeNICAM	= 1; // sec
const unsigned int kTimeLETKF	= 1; // sec
#else
const size_t kNumEnsemble       = 64;
const size_t kNumRegion         = 10;
const size_t kNumCell           = 1156;
// Assume that each lattice has 6160 KB of data (6160 = 1540 * 4)
const size_t kElementCount      = 1540;
const unsigned int kTimeNICAM   = 50; // sec
const unsigned int kTimeLETKF   = 5;  // sec
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
  MPI_Init(&argc, &argv);
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
}

void workflow(string prog_name) {
#if DEBUG
  assert(nprocs == kNumEnsemble);
  assert(universe_size >= (int)(kNumEnsemble * (1 + kNumRegion)));
#endif
  if (rank == 0) {
    cout << "PARAMETERS" << endl;
    cout << "  Universe size    : " << universe_size << endl;
    cout << "  Static proc size : " << nprocs << endl;
#ifdef SYENV_K_FX
    cout << "  Spawn gap in msec : " << spawn_gap_msec << endl;
#endif
    cout << endl;
  }

  for (size_t i = 0; i < kNumIteration; i++) {
    Time time;
    time.loop_start = gettime(MPI_COMM_WORLD);

    // run pseudo-NICAM
    spawn_nicam(prog_name, time);

    // shuffle
    shuffle(time);

    // run pseudo-LETKF
    spawn_letkf(prog_name, time);

    time.loop_finish = gettime(MPI_COMM_WORLD);
    ostringstream os1;
    os1 << "Iteration[" << i << "] ends in " << time.loop() << " sec." << endl;
    os1 << "  NICAM takes " << time.nicam() << " sec." << endl;
    os1 << "  Shuffle takes " << time.shuffle() << " sec." << endl;
    os1 << "  LETKF takes " << time.letkf() << " sec." << endl;
    print_line(os1);
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
#ifdef SYENV_K_FX
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

void spawn_letkf(string prog_name, Time& time) {
  time.letkf_start = gettime(MPI_COMM_WORLD);
  string prog_letkf("letkf");
  LETKFProc *letkfs = (LETKFProc*)calloc(letkf_task_count, sizeof(LETKFProc));
  int nmaxspawn = (universe_size - nprocs) / nprocs;
  int running = 0;

  for (size_t i = 0; i < letkf_task_count; i++) {
    if (running == nmaxspawn) {
      // Wait until a dynamic process is freed.
      for (size_t j = 0; j < i; j++) {
	if (letkfs[j].done) { continue; }
	MPI_Status stat;
	MPI_Wait(&(letkfs[j].rreq), &stat);
	free(letkfs[j].rbuf);
	disconnect(&(letkfs[j].icomm));
	letkfs[j].done = true;
	running -= 1;
      }
      assert(running < nmaxspawn);
    }

    letkfs[i].done = false;
    spawn(prog_name, prog_letkf, 1, &(letkfs[i].icomm));
    // Send data to a spawned process.
    int cc;
    int *sbuf = load_letkf_data();
    cc = MPI_Send(sbuf, (int)letkf_data_count, MPI_INT, 0, 1000,
     		  letkfs[i].icomm);
    assert(cc == MPI_SUCCESS);
    free(sbuf);

    // Receive data from a spawned process.
    letkfs[i].rbuf = allocate_data(letkf_data_count);
    cc = MPI_Irecv(letkfs[i].rbuf, (int)letkf_data_count, MPI_INT, 0, 1001,
     		   letkfs[i].icomm, &(letkfs[i].rreq));
    assert(cc == MPI_SUCCESS);

    running += 1;
  }
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

  sleep(kTimeNICAM);

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
  assert(nprocs == 1);
  MPI_Comm pcomm;
  MPI_Comm_get_parent(&pcomm);

  if (rank == 0) {
    // receive data from the spawner on rank 0
    int *buf = allocate_data(letkf_data_count);
    MPI_Status stat;
    int cc = MPI_Recv(buf, (int)letkf_data_count, MPI_INT, 0, 1000, pcomm,
    		      &stat);
    assert(cc == MPI_SUCCESS);
    free(buf);
  }

  sleep(kTimeLETKF);

  if (rank == 0) {
    int *buf = load_letkf_data();
    // send data to the spawner from rank 0
    int cc = MPI_Send(buf, (int)letkf_data_count, MPI_INT, 0, 1001, pcomm);
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
  int *buf = allocate_data(letkf_data_count);
  for (size_t i = 0; i < kNumEnsemble; i++) {
    size_t idxi = i * kElementCount;
    for (size_t k = 0; k < kElementCount; k++) {
      size_t idx = idxi + k;
      buf[idx] = (int)i;
    }
  }
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
