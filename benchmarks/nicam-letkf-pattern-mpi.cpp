/// \file
/// Benchmark program for NICAM-LETKF data access pattern
///
/// It does not use KMRNEXT but uses just MPI.
/// The number of MPI processes running this program should be
/// "kNumEnsemble x kNumRegion".
#include <iostream>
#include <sstream>
#include <cassert>
#include <ctime>
#include <mpi.h>

using namespace std;

int rank = 0;
int rank_offset = 1;

const size_t kNumIteration = 10;

const size_t kDimEnsembleData = 3;
//const size_t kDimRegionData   = 2;
//const size_t kDimCellData     = 1;

#if DEBUG
const size_t kNumEnsemble  = 2;
const size_t kNumRegion    = 10;
const size_t kNumCell      = 10;
const size_t kElementCount = 2;
#else
const size_t kNumEnsemble  = 64;
const size_t kNumRegion    = 40;
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

  double alc0_start;
  double alc0_finish;
  double alc1_start;
  double alc1_finish;

  double del0_start;
  double del0_finish;
  double del1_start;
  double del1_finish;

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
  double alc0() {
    return (alc0_finish - alc0_start) / 10E9;
  }
  double alc1() {
    return (alc1_finish - alc1_start) / 10E9;
  }
  double del0() {
    return (del0_finish - del0_start) / 10E9;
  }
  double del1() {
    return (del1_finish - del1_start) / 10E9;
  }
};

class DataStore {
public:
  DataStore(size_t ndims) : ndims_(ndims), dim_sizes_(NULL), buf_(NULL) {}

  ~DataStore() {
    if (dim_sizes_ != NULL) {
      delete[] dim_sizes_;
    }
    if (buf_ != NULL) {
      delete[] buf_;
    }
  }

  void set(const size_t *dim_sizes) {
    dim_sizes_ = new size_t[ndims_];
    for (size_t i = 0; i < ndims_; i++) {
      dim_sizes_[i] = dim_sizes[i];
    }
    ensid_ = rank / dim_sizes_[1]; // kNumRegion
    prcid_ = rank % dim_sizes_[1]; // kNumRegion

    MPI_Comm_size(MPI_COMM_WORLD, &nprocs_);
  }

  void load() {
    assert(buf_ == NULL);
    buf_size_ = dim_sizes_[2] * kElementCount; // kNumCell
    buf_ = new int[buf_size_];
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (size_t i = 0; i < buf_size_; i++) {
      buf_[i] = (int)prcid_ + 1;
    }
  }

  void dump() {
    assert(buf_ != NULL);
    for (int i = 0; i < nprocs_; i++) {
      if (rank == i) {
	for (size_t j = 0; j < buf_size_; j += kElementCount) { // kNumCell
	  ostringstream os;
	  os << "<" << ensid_ << "," << prcid_ << "," << (j/kElementCount)
	     << "> : ";
	  for (size_t k = 0; k < kElementCount; k++) {
	    os << buf_[j + k] << " ";
	  }
	  os << endl;
	  cout << os.str();
	}
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
    if (rank == 0) {
      cout << endl;
    }
  }

  void read(int **buf, size_t *buf_size) {
    assert(buf_ != NULL);
    int snd_rank = (rank + rank_offset + nprocs_) % nprocs_;
    int rcv_rank = (rank - rank_offset + nprocs_) % nprocs_;
    size_t snd_size = buf_size_;
    size_t rcv_size;
    MPI_Status st;
    MPI_Sendrecv(&snd_size, 1, MPI_LONG, snd_rank, 1000,
		 &rcv_size, 1, MPI_LONG, rcv_rank, 1000,
		 MPI_COMM_WORLD, &st);
    *buf_size = rcv_size;
    *buf = new int[*buf_size];
    MPI_Sendrecv(buf_, (int)snd_size, MPI_INT, snd_rank, 1001,
    		 *buf, (int)rcv_size, MPI_INT, rcv_rank, 1001,
    		 MPI_COMM_WORLD, &st);
  }

  void write(int *buf, size_t buf_size) {
    assert(buf_ == NULL);
    int snd_rank = (rank - rank_offset + nprocs_) % nprocs_;
    int rcv_rank = (rank + rank_offset + nprocs_) % nprocs_;
    size_t snd_size = buf_size;
    size_t rcv_size;
    MPI_Status st;
    MPI_Sendrecv(&snd_size, 1, MPI_LONG, snd_rank, 1000,
		 &rcv_size, 1, MPI_LONG, rcv_rank, 1000,
		 MPI_COMM_WORLD, &st);
    buf_size_ = rcv_size;
    buf_ = new int[buf_size_];
    MPI_Sendrecv(buf,  (int)snd_size, MPI_INT, snd_rank, 1001,
    		 buf_, (int)rcv_size, MPI_INT, rcv_rank, 1001,
    		 MPI_COMM_WORLD, &st);
  }

  size_t ensid() { return ensid_; }
  size_t prcid() { return prcid_; }

private:
  size_t ndims_;
  size_t *dim_sizes_;
  size_t ensid_;
  size_t prcid_;
  int *buf_;
  size_t buf_size_;
  int nprocs_;
};

void load_data(DataStore* ds);
void run_nicam(DataStore* inds, DataStore* outds, Time& time);
void run_letkf(DataStore* inds, DataStore* outds, Time& time);
void print_line(ostringstream& os);
double gettime();
double gettime(MPI_Comm comm);


//////////////////////////////////////////////////////////////////////////////
// Main starts from here.
//////////////////////////////////////////////////////////////////////////////
int
main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);
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
    rank_offset = nprocs / 2;
  }

  DataStore *ds0 = new DataStore(kDimEnsembleData);
  ds0->set(kEnsembleDataDimSizes);
  load_data(ds0);
#if 0
  ds0->dump();
#endif

  for (size_t i = 0; i < kNumIteration; i++) {
    Time time;
    time.loop_start = gettime();

    // run pseudo-NICAM
    time.alc0_start = gettime();
    DataStore* ds1 = new DataStore(kDimEnsembleData);
    ds1->set(kEnsembleDataDimSizes);
    time.alc0_finish = gettime();
    run_nicam(ds0, ds1, time);
    time.del0_start = gettime();
    delete ds0;
    time.del0_finish = gettime();
#if 0
    ds1->dump();
#endif

    // run pseudo-LETKF
    time.alc1_start = gettime();
    ds0 = new DataStore(kDimEnsembleData);
    ds0->set(kEnsembleDataDimSizes);
    time.alc1_finish = gettime();
    run_letkf(ds1, ds0, time);
    time.del1_start = gettime();
    delete ds1;
    time.del1_finish = gettime();
#if 0
    ds0->dump();
#endif

    time.loop_finish = gettime();
    ostringstream os1;
    os1 << "Iteration[" << i << "]," << time.loop() << endl;
    os1 << "Alc NICAM In," << time.alc0() << endl;
    os1 << "Invoke NICAM," << time.nicam_launch() << endl;
    //os1 << "NICAM,"        << time.nicam() << endl;
    os1 << "Del NICAM In," << time.del0() << endl;
    os1 << "Alc LETKF In," << time.alc1() << endl;
    os1 << "Invoke LETKF," << time.letkf_launch() << endl;
    //os1 << "LETKF,"        << time.letkf() << endl;
    os1 << "Del LETKF In," << time.del1() << endl;
    print_line(os1);
  }
  delete ds0;

  MPI_Finalize();
  return 0;
}

void load_data(DataStore* ds)
{
  ds->load();
}

void run_nicam(DataStore* inds, DataStore* outds, Time& time)
{
  time.nicam_invoke = gettime();

  MPI_Comm task_comm;
  MPI_Comm_split(MPI_COMM_WORLD, (int)inds->ensid(), rank, &task_comm);
#if DEBUG
  int nprocs_tc;
  MPI_Comm_size(task_comm, &nprocs_tc);
  assert(nprocs_tc == kNumRegion);
#endif

  time.nicam_start = gettime(task_comm);
  int *buf;
  size_t buf_size;
  inds->read(&buf, &buf_size);
#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for (size_t i = 0; i < buf_size; i++) {
    buf[i] += 1;
  }
  outds->write(buf, buf_size);
  delete[] buf;
  time.nicam_finish = gettime(task_comm);

  MPI_Comm_free(&task_comm);
  time.nicam_cleanup = gettime();
}

void letkf(int *in, int *out, size_t size, MPI_Comm comm) {
  // alltoall the input data
  int letkf_nprocs;
  int letkf_rank;
  MPI_Comm_size(comm, &letkf_nprocs);
  MPI_Comm_rank(comm, &letkf_rank);
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
	 (int)size); // kNumCell * kElementCount
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
  MPI_Alltoallv(in,     send_cnts, sdispls, MPI_INT,
		rcvbuf, recv_cnts, rdispls, MPI_INT, comm);

  // computation: just decrease the value
#ifdef _OPENMP
  #pragma omp parallel for
#endif
  for (size_t i = 0; i < (size_t)rcvbuf_siz; i++) {
    rcvbuf[i] -= 1;
  }

  // alltoall the output data
  MPI_Alltoallv(rcvbuf, recv_cnts, rdispls, MPI_INT,
		out,    send_cnts, sdispls, MPI_INT, comm);
  delete[] send_cnts;
  delete[] sdispls;
  delete[] recv_cnts;
  delete[] rdispls;
  delete[] rcvbuf;
}

void run_letkf(DataStore* inds, DataStore* outds, Time& time)
{
  time.letkf_invoke = gettime();

  MPI_Comm task_comm;
  MPI_Comm_split(MPI_COMM_WORLD, (int)inds->prcid(), rank, &task_comm);
#if DEBUG
  int nprocs_tc;
  MPI_Comm_size(task_comm, &nprocs_tc);
  assert(nprocs_tc == kNumEnsemble);
#endif

  time.letkf_start = gettime(task_comm);
  int *buf;
  size_t buf_size;
  inds->read(&buf, &buf_size);
  int *outbuf = new int[buf_size];
  letkf(buf, outbuf, buf_size, task_comm);
  delete[] buf;
  outds->write(outbuf, buf_size);
  delete[] outbuf;
  time.letkf_finish = gettime(task_comm);

  MPI_Comm_free(&task_comm);
  time.letkf_cleanup = gettime();
}

static void print_line(string& str) {
  if (kPrint && (rank == 0)) {
    cout << str << endl;
  }
}

void print_line(ostringstream& os) {
  string s = os.str();
  print_line(s);
}

double gettime() {
  return gettime(MPI_COMM_WORLD);
}

double gettime(MPI_Comm comm) {
  MPI_Barrier(comm);
  struct timespec ts;
  clock_gettime(CLOCK_REALTIME, &ts);
  return ((double) ts.tv_sec) * 10E9 + ((double) ts.tv_nsec);
}
