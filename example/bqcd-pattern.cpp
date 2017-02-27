// This program works only with KMR backend.
#include <mpi.h>
#include "kmrnext.hpp"

using namespace std;
using namespace kmrnext;

const size_t kDimGrid  = 4;
const size_t kDimWhole = 5;

// Dimension order is x, y, z and t.
// The product of all elements in kProc should be same as the number of
// MPI processes.
const size_t kGrid[kDimGrid] = {4, 4, 4, 4};
const size_t kProc[kDimGrid] = {1, 2, 1, 2};

// The product of kEnsembles and all elements in kEnsembleProc should be
// same as the product of all elements in kProc.
// Dimension order in kEnsembleProc is ensemble, x, y, z, and t.
const size_t kEnsembles = 2;
const size_t kEnsembleProc[kDimGrid] = {1, 2, 1, 1};

// It is the count of long integer that is used for data at each coordinate.
// The total amount of data of a coordinate is 'kLongCount x sizeof(long)'.
// The total amount of data of the whole grid system is the product of
// kEnsembles, all elements in kGrid, kLongCount and sizeof(long).
#if DEBUG
const size_t kLongCount = 2;
#else
const size_t kLongCount = 2 * 156;
#endif

// Test count
const size_t kIteration = 10;

void check_configuration(int nprocs, int rank);
void setup_parameters(size_t* whole_grid, size_t* coord4d, int rank);

// Data Loader class that loads data to the input DataStore.
class DataLoader : public DataStore::Loader<long> {
private:
  size_t* coord4d_;
public:
  DataLoader(size_t* coord4d) : coord4d_(coord4d) {}
  int operator()(DataStore* ds, const long& num)
  {
    size_t blocks[kDimGrid];
    for (size_t i = 0; i < kDimGrid; i++) {
      blocks[i] = kGrid[i] / kProc[i];
    }
    size_t heads[kDimGrid], tails[kDimGrid];
    for (size_t i = 0; i < kDimGrid; i++) {
      heads[i] = coord4d_[i] * blocks[i];
      tails[i]   = heads[i] + blocks[i];
    }

    // Create a data of a point
    long* buf = new long[kLongCount];
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (size_t i = 0; i < kLongCount; i++) {
      buf[i] = i+1;
    }
    Data data(static_cast<void*>(buf), sizeof(long) * kLongCount);

    // Create key and add the data
    Key key(kDimWhole);
    for (size_t i0 = heads[0]; i0 < tails[0]; i0++) {
      key.set_dim(1, i0);
      for (size_t i1 = heads[1]; i1 < tails[1]; i1++) {
	key.set_dim(2, i1);
	for (size_t i2 = heads[2]; i2 < tails[2]; i2++) {
	  key.set_dim(3, i2);
	  for (size_t i3 = heads[3]; i3 < tails[3]; i3++) {
	    key.set_dim(4, i3);
	    for (size_t e = 0; e < kEnsembles; e++) {
	      key.set_dim(0, e);

	      ds->add(key, data);
	    }
	  }
	}
      }
    }

    delete[] buf;
    return 0;
  }
};

class TaskDirac : public DataStore::Mapper {
  int world_size_;
public:
  TaskDirac(int nprocs) : world_size_(nprocs) {}

  int operator()(DataStore* inds, DataStore* outds,
		 Key& key, vector<DataPack>& dps,
		 DataStore::MapEnvironment& env)
  {
#if DEBUG
    // Check number of processes
    int nprocs;
    MPI_Comm_size(env.mpi_comm, &nprocs);
    int expect_nprocs = world_size_ / static_cast<int>(kEnsembles);
    if (nprocs != expect_nprocs) {
      cerr << "Process count should be " << expect_nprocs << "." << endl;
      KMRNext::abort(1);
    }

    // Check data count
    size_t expect_count = 1;
    size_t blocks[kDimGrid];
    for (size_t i = 0; i < kDimGrid; i++) {
      blocks[i] = kGrid[i] / kEnsembleProc[i];
      expect_count *= blocks[i];
    }
    if (dps.size() != expect_count) {
      cerr << "Data count is incorrect." << endl;
      cerr << dps.size() << " " << expect_count << endl;
      KMRNext::abort(1);
    }

    // Check region of Key
    size_t heads[kDimGrid];
    size_t tails[kDimGrid];
    for (size_t i = 0; i < kDimGrid; i++) {
      Key k = dps.at(0).key();
      heads[i] = (k.dim(i+1) / blocks[i]) * blocks[i];
      tails[i] = heads[i] + blocks[i];
    }
    for (vector<DataPack>::iterator itr = dps.begin(); itr != dps.end();
    	 itr++) {
      Key k = itr->key();
      for (size_t i = 0; i < kDimGrid; i++) {
	if (!(heads[i] <= k.dim(i+1) && k.dim(i+1) < tails[i])) {
	  cerr << "Data not expected exists." << endl;
	  KMRNext::abort(1);
	}
      }
    }
#endif
    return 0;
  }
};

class TaskMove : public DataStore::Mapper {
  int world_size_;
public:
  TaskMove(int nprocs) : world_size_(nprocs) {}

  int operator()(DataStore* inds, DataStore* outds,
		 Key& key, vector<DataPack>& dps,
		 DataStore::MapEnvironment& env)
  {
#if DEBUG
    // Check number of processes
    int nprocs;
    MPI_Comm_size(env.mpi_comm, &nprocs);
    int expect_nprocs = world_size_;
    if (nprocs != expect_nprocs) {
      cerr << "Process count should be " << expect_nprocs << "." << endl;
      KMRNext::abort(1);
    }

    // Check data count
    size_t expect_count = kEnsembles;
    size_t blocks[kDimGrid];
    for (size_t i = 0; i < kDimGrid; i++) {
      blocks[i] = kGrid[i] / kProc[i];
      expect_count *= blocks[i];
    }
    if (dps.size() != expect_count) {
      cerr << "Data count is incorrect." << endl;
      cerr << dps.size() << " " << expect_count << endl;
      KMRNext::abort(1);
    }

    // Check region of Key
    size_t heads[kDimGrid];
    size_t tails[kDimGrid];
    for (size_t i = 0; i < kDimGrid; i++) {
      Key k = dps.at(0).key();
      heads[i] = (k.dim(i+1) / blocks[i]) * blocks[i];
      tails[i] = heads[i] + blocks[i];
    }
    for (vector<DataPack>::iterator itr = dps.begin(); itr != dps.end();
    	 itr++) {
      Key k = itr->key();
      for (size_t i = 0; i < kDimGrid; i++) {
	if (!(heads[i] <= k.dim(i+1) && k.dim(i+1) < tails[i])) {
	  cerr << "Data not expected exists." << endl;
	  KMRNext::abort(1);
	}
      }
    }
#endif
    return 0;
  }
};


//////////////////////////////////////////////////////////////////////////////
// Main starts from here.
//////////////////////////////////////////////////////////////////////////////
int
main(int argc, char **argv)
{
  KMRNext *next = KMRNext::init(argc, argv);
  //next->enable_profile();

  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  check_configuration(nprocs, rank);

  // Prepare execution
  size_t whole_grid[kDimWhole];
  size_t coord4d[kDimGrid];
  setup_parameters(whole_grid, coord4d, rank);

  DataLoader loader(coord4d);
  View dView(kDimWhole);
  View dSplit(kDimWhole);
  {
    long view_flags[kDimWhole];  // A, N, N, N, N
    long spli_flags[kDimWhole];  // A, n, n, n, n
    view_flags[0] = View::SplitAll;
    spli_flags[0] = View::SplitAll;
    for (size_t i = 1; i < kDimWhole; i++) {
      view_flags[i] = View::SplitNone;
      spli_flags[i] = kEnsembleProc[i-1];
    }
    dView.set(view_flags);
    dSplit.set(spli_flags);
  }

  TaskMove task_move(nprocs);
  View mView(kDimWhole);
  View mSplit(kDimWhole);
  {
    long view_flags[kDimWhole];  // N, N, N, N, N
    long spli_flags[kDimWhole];  // A, n, n, n, n
    view_flags[0] = View::SplitNone;
    spli_flags[0] = View::SplitNone;
    for (size_t i = 1; i < kDimWhole; i++) {
      view_flags[i] = View::SplitNone;
      spli_flags[i] = kProc[i-1];
    }
    mView.set(view_flags);
    mSplit.set(spli_flags);
  }

  // Create a DS and add data to it on each process
  DataStore* ds0 = next->create_ds(kDimWhole);
  ds0->set(whole_grid);
  ds0->load_parallel(loader);

  for (size_t i = 0; i < kIteration; i++) {
    // Run dirac solver
    TaskDirac task_dirac(nprocs);
    ds0->set_split(dSplit);
    ds0->map(task_dirac, dView);
#if DEBUG
    if (ds0->collated()) {
      if (rank == 0) {
	cerr << "collated at task_dirac." << endl;
      }
    }
#endif

    // Restore to the original data layout
    ds0->set_split(mSplit);
    ds0->map(task_move, mView);
#if DEBUG
    if (ds0->collated()) {
      if (rank == 0) {
	cerr << "collated at task_move." << endl;
      }
    }
#endif
  }

  delete ds0;
  KMRNext::finalize();
  return 0;
}

void check_configuration(int nprocs, int rank)
{
  int expect_nprocs = 1;
  for (size_t i = 0; i < kDimGrid; i++) {
    expect_nprocs *= static_cast<int>(kProc[i]);
  }
  if (nprocs != expect_nprocs) {
    if (rank == 0) {
      cerr << "The number of processes is illegal." << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    KMRNext::abort(1);
  }

  for (size_t i = 0; i < kDimGrid; i++) {
    if (kGrid[i] % kProc[i] != 0) {
      if (rank == 0) {
	cerr << "Division of " << i << "'th dimension is illegal." << endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
      KMRNext::abort(1);
    }
  }

  expect_nprocs = kEnsembles;
  for (size_t i = 0; i < kDimGrid; i++) {
    expect_nprocs *= static_cast<int>(kEnsembleProc[i]);
  }
  if (nprocs != expect_nprocs) {
    if (rank == 0) {
      cerr << "The number of processes in the ensemble computation is illegal."
	   << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    KMRNext::abort(1);
  }
}

void setup_parameters(size_t* whole_grid, size_t* coord4d, int rank)
{
  whole_grid[0] = kEnsembles;
  for (size_t i = 1; i < kDimWhole; i++) {
    whole_grid[i] = kGrid[i-1];
  }

  size_t tmp = rank;
  for (size_t i = 0; i < kDimGrid; i++) {
    coord4d[i] = tmp % kProc[i];
    tmp /= kProc[i];
  }
}
