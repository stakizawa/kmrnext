#include <iostream>
#include "kmrnext.hpp"
#ifdef BACKEND_KMR
#include <mpi.h>
#endif
#include <cassert>
/// This tests creating multiple instances of KMRNext.

using namespace std;

int rank = 0;
int nprocs = 1;

const bool kPrint = true;

const size_t kN = 4;
const size_t kDim1 = 1;
const size_t kDim2 = 2;

void session(kmrnext::KMRNext *next, long *vals, size_t nvals);
void print_data(long *vals, size_t nvals);


//////////////////////////////////////////////////////////////////////////////
// Main starts from here.
//////////////////////////////////////////////////////////////////////////////
int
main(int argc, char **argv) {
#ifdef BACKEND_KMR
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
#endif

  long *vals = new long[kN];
  for (size_t i = 0; i < kN; i++) {
    vals[i] = rank;
  }
  if (kPrint) {
    print_data(vals, kN);
  }

  ///////////  The first KMRNext
  kmrnext::KMRNext *next0 = kmrnext::KMRNext::init(argc, argv);
  session(next0, vals, kN);
  kmrnext::KMRNext::finalize();

  ///////////  The second KMRNext
  kmrnext::KMRNext *next1 = kmrnext::KMRNext::init();
  session(next1, vals, kN);
  kmrnext::KMRNext::finalize();

  delete[] vals;

#ifdef BACKEND_KMR
  MPI_Finalize();
#endif
  return 0;
}


class IntLoader2 : public kmrnext::DataStore::Loader<long> {
private:
  long *vals_;
  size_t nvals_;
public:
  IntLoader2(long *vals, size_t nvals) : vals_(vals), nvals_(nvals) {}

  int operator()(kmrnext::DataStore* ds, const long& i) {
    kmrnext::Key key(kDim2);
    key.set_dim(0, i);
    for (size_t j = 0; j < nvals_; j++) {
      key.set_dim(1, j);
      kmrnext::Data d(&vals_[j], sizeof(long));
      ds->add(key, d);
    }
    return 0;
  }
};

class Summarizer : public kmrnext::DataStore::Mapper {
private:
  long ndata_;
public:
  Summarizer(long ndata) : ndata_(ndata) {}
  int operator()(kmrnext::DataStore* inds, kmrnext::DataStore* outds,
		 kmrnext::Key& key, vector<kmrnext::DataPack>& dps,
		 kmrnext::DataStore::MapEnvironment& env) {
    assert(dps.size() == static_cast<size_t>(ndata_));
    long sum = 0;
    for (size_t i = 0; i < dps.size(); i++) {
      kmrnext::DataPack& dp = dps.at(i);
      long v = *static_cast<long*>(dp.data().value());
      sum += v;
    }
    kmrnext::Data d(&sum, sizeof(long));
    outds->add(key, d);
    return 0;
  }
};

class Copier : public kmrnext::DataStore::Mapper {
private:
  long *ary_;
  size_t nary_;
public:
  Copier(long *ary, size_t nary) : ary_(ary), nary_(nary) {}
  int operator()(kmrnext::DataStore* inds, kmrnext::DataStore* outds,
		 kmrnext::Key& key, vector<kmrnext::DataPack>& dps,
		 kmrnext::DataStore::MapEnvironment& env) {
    assert(dps.size() == nary_);
    for (size_t i = 0; i < nary_; i++) {
      kmrnext::DataPack& dp = dps.at(i);
      ary_[i] = *static_cast<long*>(dp.data().value());
    }
    return 0;
  }
};

void session(kmrnext::KMRNext *next, long *vals, size_t nvals) {
  // Load the above array to DataStore
  kmrnext::DataStore *ds0 = next->create_ds(kDim2);
  {
    size_t size2[kDim2] = {static_cast<size_t>(nprocs), kN};
    ds0->set(size2);
    vector<long> ints;
    for (int i = 0; i < nprocs; i++) {
      ints.push_back(i);
    }
    IntLoader2 loader(vals, kN);
    ds0->load_integers(ints, loader);
  }
  // Apply map to the DataStore
  kmrnext::DataStore *ds1 = next->create_ds(kDim1);
  {
    size_t size1[kDim1] = {kN};
    ds1->set(size1);
    Summarizer mapper(nprocs);
#ifdef BACKEND_KMR
    kmrnext::View split(kDim2);
    bool split_flag[kDim2] = {false, true};
    split.set(split_flag);
    ds0->set_split(split);
#endif
    kmrnext::View view(kDim2);
    bool view_flag[kDim2] = {false, true};
    view.set(view_flag);
    ds0->map(mapper, view, ds1);
  }
  delete ds0;

  // Apply map to copy data
  long val;
  {
    Copier mapper(&val, 1);
#ifdef BACKEND_KMR
    kmrnext::View split(kDim1);
    bool split_flag[kDim1] = {true};
    split.set(split_flag);
    ds1->set_split(split);
#endif
    kmrnext::View view(kDim1);
    bool view_flag[kDim1] = {true};
    view.set(view_flag);
    ds1->map(mapper, view);
  }
  if (kPrint) {
    print_data(&val, 1);
  }
  delete ds1;
}

void print_data(long *vals, size_t nvals) {
#ifdef BACKEND_KMR
  for (int j = 0; j < static_cast<int>(nprocs); j++) {
    if (rank == j) {
      for (size_t i = 0; i < nvals; i++) {
	cout << vals[i] << " ";
      }
      cout << endl;
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
#else
  for (size_t i = 0; i < nvals; i++) {
    cout << vals[i] << " ";
  }
  cout << endl;
#endif
}
