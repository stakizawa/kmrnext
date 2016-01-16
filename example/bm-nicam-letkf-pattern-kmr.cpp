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
const unsigned int kTimeNICAM	= 0; // sec
const unsigned int kTimeLETKF	= 0; // sec
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

int rank		= 0;
int nprocs		= 1;
size_t nicam_data_count = 0;
size_t letkf_data_count = 0;

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

  nicam_data_count = kNumRegion * kNumCell * kElementCount;
  int div_nicam = (int)kNumEnsemble / nprocs;
  int rem_nicam = (int)kNumEnsemble % nprocs;
  nicam_data_count *= (rank < rem_nicam)? div_nicam + 1 : div_nicam;

  letkf_data_count = kNumEnsemble * kNumCell * kElementCount;
  int div_letkf = (int)kNumRegion / nprocs;
  int rem_letkf = (int)kNumRegion % nprocs;
  letkf_data_count *= (rank < rem_letkf)? div_letkf + 1 : div_letkf;
}

void workflow(string prog_name) {
#if DEBUG
  assert(nprocs == kNumEnsemble);
  int *usizep;
  int uflag;
  MPI_Comm_get_attr(MPI_COMM_WORLD, MPI_UNIVERSE_SIZE, &usizep, &uflag);
  int univ_size = *usizep;
  assert(univ_size >= (int)(kNumEnsemble * (1 + kNumRegion)));
#endif

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
}

void spawn_nicam(string prog_name, Time& time) {
  time.nicam_start = gettime(MPI_COMM_WORLD);
  string prog_nicam("nicam");
  MPI_Comm icomm;
  spawn(prog_name, prog_nicam, kNumRegion, &icomm);

  int cc;
  int *sbuf = load_nicam_data();
  cc = MPI_Send(sbuf, (int)nicam_data_count, MPI_INT, 0, 1000, icomm);
  assert(cc == MPI_SUCCESS);

  int *rbuf = allocate_data(nicam_data_count);
  MPI_Status stat;
  cc = MPI_Recv(rbuf, (int)nicam_data_count, MPI_INT, 0, 1001, icomm, &stat);
  assert(cc == MPI_SUCCESS);

  MPI_Comm_disconnect(&icomm);
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
  time.shuffle_finish = gettime(MPI_COMM_WORLD);
}

void spawn_letkf(string prog_name, Time& time) {
  time.letkf_start = gettime(MPI_COMM_WORLD);
  string prog_letkf("letkf");
  MPI_Comm icomm;
  spawn(prog_name, prog_letkf, kNumEnsemble, &icomm);

  int cc;
  int *sbuf = load_letkf_data();
  cc = MPI_Send(sbuf, (int)letkf_data_count, MPI_INT, 0, 1000, icomm);
  assert(cc == MPI_SUCCESS);

  int *rbuf = allocate_data(letkf_data_count);
  MPI_Status stat;
  cc = MPI_Recv(rbuf, (int)letkf_data_count, MPI_INT, 0, 1001, icomm, &stat);
  assert(cc == MPI_SUCCESS);

  MPI_Comm_disconnect(&icomm);
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
  assert(nprocs == kNumEnsemble);
  MPI_Comm pcomm;
  MPI_Comm_get_parent(&pcomm);

  if (rank == 0) {
    // receive data from the spawner on rank 0
    int *buf = allocate_data(letkf_data_count);
    MPI_Status stat;
    int cc = MPI_Recv(buf, (int)letkf_data_count, MPI_INT, 0, 1000, pcomm,
		      &stat);
    assert(cc == MPI_SUCCESS);
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
  size_t ensembles =
    nicam_data_count / (kNumRegion * kNumCell * kElementCount);
  for (size_t h = 0; h < ensembles; h++) {
    size_t idxh = h * kNumRegion * kNumCell * kElementCount;
    for (size_t i = 0; i < kNumRegion; i++) {
      size_t idxi = i * kNumCell * kElementCount;
      for (size_t j = 0; j < kNumCell; j++) {
	size_t idxj = j * kElementCount;
	for (size_t k = 0; k < kElementCount; k++) {
	  size_t idx = idxh + idxi + idxj + k;
	  buf[idx] = (int)i;
	}
      }
    }
  }
  return buf;
}

int* load_letkf_data() {
  int *buf = allocate_data(letkf_data_count);
  size_t regions =
    letkf_data_count / (kNumEnsemble * kNumCell * kElementCount);
  for (size_t h = 0; h < regions; h++) {
    size_t idxh = h * kNumEnsemble * kNumCell * kElementCount;
    for (size_t i = 0; i < kNumEnsemble; i++) {
      size_t idxi = i * kNumCell * kElementCount;
      for (size_t j = 0; j < kNumCell; j++) {
	size_t idxj = j * kElementCount;
	for (size_t k = 0; k < kElementCount; k++) {
	  size_t idx = idxh + idxi + idxj + k;
	  buf[idx] = (int)i;
	}
      }
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
