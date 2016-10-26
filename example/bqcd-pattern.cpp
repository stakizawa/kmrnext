#include <iostream>
#include <mpi.h>
#include "kmrnext.hpp"

using namespace std;
using namespace kmrnext;

const size_t kDimSize = 4;

const size_t kGrid[kDimSize] = {4, 4, 4, 4};
const size_t kProc[kDimSize] = {1, 2, 1, 2};

// The product of kEnsembles and all elements in kEnsembleProc should be
// same as the product of all elements in kProc
const size_t kEnsembles = 2;
const size_t kEnsembleProc[kDimSize] = {1, 2, 1, 1};

void check_configuration(int nprocs, int rank);
void zeroize_ds(DataStore* ds);
void calc_region(const int rank, const size_t siz, size_t* head, size_t* tail,
		 const size_t *grid, const size_t *split);
void check_data_location_4d(DataStore* ds);
void check_data_location_5d(DataStore* ds);
void print_ds(DataStore* ds);

class LoadDataMapper : public DataStore::Mapper {
  size_t* heads_;
  size_t* tails_;
  size_t  size_;
  long    value_;
public:
  LoadDataMapper(size_t* heads, size_t* tails, size_t siz, long data)
    : heads_(heads), tails_(tails), size_(siz), value_(data) {}

  int operator()(DataStore* inds, DataStore* outds,
		 Key& key, vector<DataPack>& dps,
		 DataStore::MapEnvironment& env)
  {
    // Check data count
    size_t expected = 1;
    for (size_t i = 0; i < size_; i++) {
      expected *= tails_[i] - heads_[i];
    }
    if (expected != dps.size()) {
      cerr << "Data count is not correct." << endl;
      KMRNext::abort(1);
    }

    // Check number of processes
    int nprocs;
    MPI_Comm_size(env.mpi_comm, &nprocs);
    if (nprocs != 1) {
      cerr << "Process count should be 1." << endl;
      KMRNext::abort(1);
    }

    Data data(&value_, sizeof(long));
    for (vector<DataPack>::iterator itr = dps.begin(); itr != dps.end();
	 itr++) {
      // Check key range
      for (size_t i = 0; i < size_; i++) {
	size_t key_idx = itr->key().dim(i);
	if (!(heads_[i] <= key_idx && key_idx < tails_[i])) {
	  cerr << "Key is out of range." << endl;
	  KMRNext::abort(1);
	}
      }
      outds->add(itr->key(), data);
    }

    return 0;
  }
};

class CopyMapper : public DataStore::Mapper {
  int world_size_;
public:
  CopyMapper(int nprocs) : world_size_(nprocs) {}

  int operator()(DataStore* inds, DataStore* outds,
		 Key& key, vector<DataPack>& dps,
		 DataStore::MapEnvironment& env)
  {
    // Check number of processes
    int nprocs;
    MPI_Comm_size(env.mpi_comm, &nprocs);
    if (nprocs != world_size_) {
      cerr << "Process count should be 1." << endl;
      KMRNext::abort(1);
    }

    for (vector<DataPack>::iterator itr = dps.begin(); itr != dps.end();
    	 itr++) {
      Key nkey(inds->size() + 1);
      for (size_t i = 0; i < itr->key().size(); i++) {
	nkey.set_dim(i + 1, itr->key().dim(i));
      }
      for (size_t i = 0; i < kEnsembles; i++) {
	nkey.set_dim(0, i);
	Data dat(itr->data().value(), itr->data().size());
	outds->add(nkey, dat);
      }
    }

    return 0;
  }
};

class DiracMapper : public DataStore::Mapper {
  int world_size_;
public:
  DiracMapper(int nprocs) : world_size_(nprocs) {}

  int operator()(DataStore* inds, DataStore* outds,
		 Key& key, vector<DataPack>& dps,
		 DataStore::MapEnvironment& env)
  {
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
    size_t blocks[kDimSize];
    for (size_t i = 0; i < kDimSize; i++) {
      blocks[i] = kGrid[i] / kEnsembleProc[i];
      expect_count *= blocks[i];
    }
    if (dps.size() != expect_count) {
      cerr << "Data count is incorrect." << endl;
      KMRNext::abort(1);
    }

    // Check region of Key
    size_t heads[kDimSize];
    size_t tails[kDimSize];
    for (size_t i = 0; i < kDimSize; i++) {
      Key k = dps.at(0).key();
      heads[i] = (k.dim(i+1) / blocks[i]) * blocks[i];
      tails[i] = heads[i] + blocks[i];
    }
    for (vector<DataPack>::iterator itr = dps.begin(); itr != dps.end();
    	 itr++) {
      Key k = itr->key();
      for (size_t i = 0; i < kDimSize; i++) {
	if (!(heads[i] <= k.dim(i+1) && k.dim(i+1) < tails[i])) {
	  cerr << "Data not expected exists." << endl;
	  KMRNext::abort(1);
	}
      }
    }

    // Set data
    long value = dps.at(0).key().dim(0);
    Data dat(&value, sizeof(long));
    for (vector<DataPack>::iterator itr = dps.begin(); itr != dps.end();
    	 itr++) {
      outds->add(itr->key(), dat);
    }

    return 0;
  }
};

class MoveMapper : public DataStore::Mapper {
  int world_size_;
public:
  MoveMapper(int nprocs) : world_size_(nprocs) {}

  int operator()(DataStore* inds, DataStore* outds,
		 Key& key, vector<DataPack>& dps,
		 DataStore::MapEnvironment& env)
  {
    // Check number of processes
    int nprocs;
    MPI_Comm_size(env.mpi_comm, &nprocs);
    if (nprocs != world_size_) {
      cerr << "Process count should be " << world_size_ << "." << endl;
      KMRNext::abort(1);
    }

    // Check data count
    size_t expect_count = kEnsembles;
    size_t blocks[kDimSize];
    for (size_t i = 0; i < kDimSize; i++) {
      blocks[i] = kGrid[i] / kProc[i];
      expect_count *= blocks[i];
    }
    if (dps.size() != expect_count) {
      cerr << "Data count is incorrect." << endl;
      KMRNext::abort(1);
    }

    // Check region of Key
    size_t heads[kDimSize];
    size_t tails[kDimSize];
    for (size_t i = 0; i < kDimSize; i++) {
      Key k = dps.at(0).key();
      heads[i] = (k.dim(i+1) / blocks[i]) * blocks[i];
      tails[i] = heads[i] + blocks[i];
    }
    for (vector<DataPack>::iterator itr = dps.begin(); itr != dps.end();
    	 itr++) {
      Key k = itr->key();
      for (size_t i = 0; i < kDimSize; i++) {
	if (!(heads[i] <= k.dim(i+1) && k.dim(i+1) < tails[i])) {
	  cerr << "Data not expected exists." << endl;
	  KMRNext::abort(1);
	}
      }
    }

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

  int nprocs;
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  check_configuration(nprocs, rank);

  // Calculate the region in the gird this process takes care
  // TODO move renge calculation to LoadMapper
  size_t *reg_heads = new size_t[kDimSize];
  size_t *reg_tails = new size_t[kDimSize];
  calc_region(rank, kDimSize, reg_heads, reg_tails, kGrid, kProc);
#if 0  // for debug
  if (rank == 0) {
    for (int i = 0; i < nprocs; i++) {
      calc_region(i, kDimSize, reg_heads, reg_tails, kGrid, kProc);
      cout << i << ": ";
      for (size_t j = 0; j < kDimSize; j++) {
	cout << "(" << reg_heads[j] << ":" << reg_tails[j] << ") ";
      }
      cout << endl;
    }
  }
#endif

  // Create a DS and set data on each process
  DataStore* ds0 = next->create_ds(kDimSize);
  ds0->set(kGrid);
  zeroize_ds(ds0);
  View load_split = View(kDimSize);
  View load_view = View(kDimSize);
  {
    load_split.set(reinterpret_cast<long*>(const_cast<size_t*>(kProc)));
    load_view.set(reinterpret_cast<long*>(const_cast<size_t*>(kProc)));
  }
#if 0 // TODO use copy constructor
  View load_view = load_split;
#endif
  LoadDataMapper ldm(reg_heads, reg_tails, kDimSize, rank);
  ds0->set_split(load_split);
  ds0->map(ldm, load_view);
  delete[] reg_heads;
  delete[] reg_tails;
  check_data_location_4d(ds0);

  // Create a DS for ensemble
  DataStore* ds1 = next->create_ds(kDimSize + 1);
  size_t* ds1_ary = new size_t[kDimSize + 1];
  ds1_ary[0] = kEnsembles;
  for (size_t i = 0; i < kDimSize; i++) {
    ds1_ary[i+1] = kGrid[i];
  }
  ds1->set(ds1_ary);

  // Copy data to ensemble DS
  CopyMapper copy(nprocs);
  View copy_view = View(kDimSize);
  {
    for (size_t i = 0; i < kDimSize; i++) {
      copy_view.set_dim(i, View::SplitNone);
    }
  }
  ds0->map(copy, copy_view, ds1);
  check_data_location_5d(ds1);

  // Run Dirac caluculation
  DiracMapper dirac(nprocs);
  View dirac_split(kDimSize + 1);
  View dirac_view(kDimSize + 1);
  {
    dirac_split.set_dim(0, View::SplitAll);
    dirac_view.set_dim(0, View::SplitAll);
    for (size_t i = 1; i < dirac_split.size(); i++) {
      dirac_split.set_dim(i, kEnsembleProc[i-1]);
      dirac_view.set_dim(i, View::SplitNone);
    }
  }
  ds1->set_split(dirac_split);
  ds1->map(dirac, dirac_view);

  // Restore to the original data layout
  MoveMapper move(nprocs);
  View move_split(kDimSize + 1);
  View move_view(kDimSize + 1);
  {
    move_split.set_dim(0, View::SplitNone);
    for (size_t i = 0; i < kDimSize; i++) {
      move_split.set_dim(i+1, kProc[i]);
    }
    for (size_t i = 0; i < kDimSize + 1; i++) {
      move_view.set_dim(i, View::SplitNone);
    }
  }
  ds1->set_split(move_split);
  ds1->map(move, move_view);

  KMRNext::finalize();
  return 0;
}


void check_configuration(int nprocs, int rank)
{
  int expect_nprocs = 1;
  for (size_t i = 0; i < kDimSize; i++) {
    expect_nprocs *= static_cast<int>(kProc[i]);
  }
  if (nprocs != expect_nprocs) {
    if (rank == 0) {
      cerr << "The number of processes is illegal." << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
    KMRNext::abort(1);
  }

  for (size_t i = 0; i < kDimSize; i++) {
    if (kGrid[i] % kProc[i] != 0) {
      if (rank == 0) {
	cerr << "The number of processes is illegal." << endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
      KMRNext::abort(1);
    }
  }

  size_t ens_div = kEnsembles;
  while (ens_div != 1) {
    if (ens_div % 2 != 0) {
      if (rank == 0) {
	cerr << "Number of ensembles should be power of 2." << endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
      KMRNext::abort(1);
    }
    ens_div /= 2;
  }

  size_t prod_proc = 1;
  size_t prod_ens_proc = kEnsembles;
  for (size_t i = 0; i < kDimSize; i++) {
    prod_proc *= kProc[i];
    prod_ens_proc *= kEnsembleProc[i];
  }
  if (prod_proc != prod_ens_proc) {
      if (rank == 0) {
	cerr << "The product of kEnsemble and elements in kEnsembleProc "
	     << "should be same as the product of elements in kProc"<< endl;
      }
      MPI_Barrier(MPI_COMM_WORLD);
      KMRNext::abort(1);
  }
}

class ZeroLoader : public DataStore::Loader<long> {
public:
  int operator()(DataStore* ds, const long& num)
  {
    Key key(ds->size());
    set_data(ds, key, 0, ds->size());
    return 0;
  }
private:
  void set_data(DataStore* ds, Key& key, size_t cur_depth, size_t max_depth)
  {
    if (cur_depth < max_depth) {
      for (size_t i = 0; i < ds->dim(cur_depth); i++) {
	key.set_dim(cur_depth, i);
	set_data(ds, key, cur_depth + 1, max_depth);
      }
    } else {
      long val = 0;
      Data data(&val, sizeof(long));
      ds->add(key, data);
    }
  }
};

void zeroize_ds(DataStore* ds)
{
  ZeroLoader zl;
  vector<long> dlst;
  dlst.push_back(0);
  ds->load_integers(dlst, zl);
}

void calc_region(const int rank, const size_t siz, size_t* head, size_t* tail,
		 const size_t *grid, const size_t *split)
{
  size_t rank_s = static_cast<size_t>(rank);
  // Order: x, y, z, t
  for (size_t i = 0; i < siz; i++) {
    size_t prod_prev = 1;
    for (long j = i - 1; j >= 0; j--) {
      prod_prev *= split[j];
    }
    head[i] = (rank_s / prod_prev) % split[i];
  }

  size_t *blocks = new size_t[siz];
  for (size_t i = 0; i < siz; i++) {
    blocks[i] = grid[i] / split[i];
  }

  for (size_t i = 0; i < siz; i++) {
    head[i] *= blocks[i];
    tail[i] = head[i] + blocks[i];
  }
  delete[] blocks;
}

int calc_owner(const size_t *idx, const size_t *grid, size_t siz)
{
  size_t owner = 0;
  for (size_t i = 0; i < siz; i++) {
    size_t prod_prev = 1;
    for (long j = i - 1; j >= 0; j--) {
      prod_prev *= grid[j];
    }
    owner += (idx[i] * prod_prev);
  }
  return static_cast<int>(owner);
}

void check_data_location_4d(DataStore* ds)
{
  size_t blocks[kDimSize];
  for (size_t i = 0; i < 4; i++) {
    blocks[i] = kGrid[i] / kProc[i];
  }

  for (size_t i = 0; i < ds->dim(0); i++) {
    for (size_t j = 0; j < ds->dim(1); j++) {
      for (size_t k = 0; k < ds->dim(2); k++) {
	for (size_t l = 0; l < ds->dim(3); l++) {
	  // To cache data, get() is required.
	  Key key(4);
	  size_t ary0[4] = {i, j, k, l};
	  key.set(ary0);
	  DataPack dp = ds->get(key);

	  // Calculate owner
	  size_t ary1[4];
	  for (size_t m = 0; m < 4; m++) {
	    ary1[m] = ary0[m] / blocks[m];
	  }
	  int owner = calc_owner(ary1, kProc, 4);

	  DataElement *de = ds->data_element_at(key);
	  if (de->owner() != owner) {
	    cerr << "Owner is not correct." << endl;
	    KMRNext::abort(1);
	  }
	}
      }
    }
  }
}

void check_data_location_5d(DataStore* ds)
{
  size_t blocks[kDimSize];
  for (size_t i = 0; i < 4; i++) {
    blocks[i] = kGrid[i] / kProc[i];
  }

  for (size_t i = 0; i < ds->dim(1); i++) {
    for (size_t j = 0; j < ds->dim(2); j++) {
      for (size_t k = 0; k < ds->dim(3); k++) {
	for (size_t l = 0; l < ds->dim(4); l++) {

	  // Ensemble dimension
	  for (size_t m = 0; m < kEnsembles; m++) {
	    // To cache data, get() is required.
	    Key key(5);
	    size_t ary0[5] = {m, i, j, k, l};
	    key.set(ary0);
	    DataPack dp = ds->get(key);

	    // Calculate owner
	    size_t ary1[4];
	    for (size_t n = 0; n < 4; n++) {
	      ary1[n] = ary0[n+1] / blocks[n];
	    }
	    int owner = calc_owner(ary1, kProc, 4);

	    DataElement *de = ds->data_element_at(key);
	    if (de->owner() != owner) {
	      cerr << "Owner is not correct." << endl;
	      KMRNext::abort(1);
	    }
	  }

	}
      }
    }
  }
}

class DataPrinter : public DataPack::Dumper {
public:
  string operator()(DataPack& dp)
  {
    ostringstream os;
    os << dp.key().to_string() << " : ";
    long *data = static_cast<long*>(dp.data().value());
    os << *data << endl;
    return os.str();
  }
};

void print_ds(DataStore* ds) {
  DataPrinter dp;
  string str = ds->dump(dp);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    cout << str << endl;
  }
}
