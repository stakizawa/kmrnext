// This test tests DataStore class for KMR backend.
#include <gtest/gtest.h>
#include <sstream>
#include <cassert>
#include <mpi.h>
#include "kmrnext.hpp"

extern kmrnext::KMRNext *gNext;

namespace {

  class DataLoader1D : public kmrnext::DataStore::Loader<int> {
    size_t size_;
  public:
    DataLoader1D(size_t siz) : size_(siz) {}

    int operator()(kmrnext::DataStore *ds, const int& num)
    {
      kmrnext::Data data((void*)&num, sizeof(int));
      for (size_t i = 0; i < size_; i++) {
	kmrnext::Key key = ds->index_to_key(i);
	ds->add(key, data);
      }
      return 0;
    }
  };

  class KMRDataStoreTest : public ::testing::Test {
  protected:
    KMRDataStoreTest() {
      MPI_Comm_rank(gNext->kmr()->comm, &rank);
      MPI_Comm_size(gNext->kmr()->comm, &nprocs);

      ds2_array_ = new size_t[2];
      ds2_array_[0] = 4;
      ds2_array_[1] = 4;
      std::vector<int> ds2_vec;
      for (size_t i = 0; i < ds2_array_[0]; i++) {
	ds2_vec.push_back(i + 1);
      }
      ds2_owners_ = new int[4];
      init_owners(ds2_owners_, 4);
      DataLoader1D ds2_loader1d(ds2_array_[1]);
      ds2_ = new kmrnext::DataStore(2, gNext);
      ds2_->set(ds2_array_);
      ds2_->load_array(ds2_vec, ds2_loader1d);

      size_t ary_k2_00[2] = {0,0};
      init_key(&k2_00_, 2, ary_k2_00);
      size_t ary_k2_11[2] = {1,1};
      init_key(&k2_11_, 2, ary_k2_11);
      size_t ary_k2_22[2] = {2,2};
      init_key(&k2_22_, 2, ary_k2_22);
      size_t ary_k2_33[2] = {3,3};
      init_key(&k2_33_, 2, ary_k2_33);
    }

    void init_key(kmrnext::Key **kp, size_t size, size_t *array) {
      *kp = new kmrnext::Key(size);
      (*kp)->set(array);
    }

    void init_owners(int *owners, size_t size) {
      size_t quotient = size / nprocs;
      size_t remain   = size % nprocs;
      for (size_t i = 0; i < (size_t)nprocs; i++) {
	size_t start = i * quotient + ((i < remain)? i : remain);
	size_t end = start + quotient + ((i < remain)? 1 : 0);
	if (start >= size) { break; }
	for (size_t j = start; j < end; j++) {
	  owners[j] = i;
	}
      }
    }

    virtual ~KMRDataStoreTest() {
      delete ds2_array_;
      delete ds2_;
      delete ds2_owners_;
      delete k2_00_;
      delete k2_11_;
      delete k2_22_;
      delete k2_33_;
    }

    int rank;
    int nprocs;

    size_t *ds2_array_;       // {4,4}
    kmrnext::DataStore *ds2_; // Dimension: {4,4}
                              // Data: <0,x> => 1, <1,x> => 2, <2,x> => 3
                              //       <3,x> => 4
    int *ds2_owners_;         // size is 4, calculated at runtime

    kmrnext::Key *k2_00_;     // <0,0>
    kmrnext::Key *k2_11_;     // <0,0>
    kmrnext::Key *k2_22_;     // <2,2>
    kmrnext::Key *k2_33_;     // <3,3>
  };

  TEST_F(KMRDataStoreTest, Load_array) {
    // also tests ds.get()
    size_t ds_size[2] = {4,4};
    int owners[4];
    DataLoader1D loader1d(ds_size[1]);

    // Correctly load data by each process
    std::vector<int> vec0;
    for (size_t i = 0; i < ds_size[0]; i++) {
      vec0.push_back(i + 1);
    }
    init_owners(owners, 4);
    kmrnext::DataStore ds0(2, gNext);
    ds0.set(ds_size);
    ds0.load_array(vec0, loader1d);
    EXPECT_EQ(1, *(int*)ds0.get(*k2_00_).data()->value());
    EXPECT_EQ(2, *(int*)ds0.get(*k2_11_).data()->value());
    EXPECT_EQ(3, *(int*)ds0.get(*k2_22_).data()->value());
    EXPECT_EQ(4, *(int*)ds0.get(*k2_33_).data()->value());
    EXPECT_EQ(owners[0], ds0.get(*k2_00_).data()->owner());
    EXPECT_EQ(owners[1], ds0.get(*k2_11_).data()->owner());
    EXPECT_EQ(owners[2], ds0.get(*k2_22_).data()->owner());
    EXPECT_EQ(owners[3], ds0.get(*k2_33_).data()->owner());
  }

  TEST_F(KMRDataStoreTest, Get_view) {
    // Check owners when View<T, F> is given.
    kmrnext::View v0(2);
    bool flags0[2] = {true, false};
    v0.set(flags0);
    std::vector<kmrnext::DataPack> *vec0 = ds2_->get(v0, *k2_00_);
    EXPECT_EQ(4, (int)vec0->size());
    for (std::vector<kmrnext::DataPack>:: iterator itr = vec0->begin();
	 itr != vec0->end(); itr++) {
      EXPECT_EQ(ds2_owners_[0], (*itr).data()->owner());
    }
    delete vec0;
  }

#if 0
  TEST_F(KMRDataStoreTest, Set_from) {
  }
#endif

#if 0
  TEST_F(KMRDataStoreTest, Split_to) {
  }
#endif


  // A mapper class that just copies a data in the input DataStore to the
  // output DataStore.
  class Copier0 : public kmrnext::DataStore::Mapper {
  public:
    int operator()(kmrnext::DataStore *inds, kmrnext::DataStore *outds,
		   kmrnext::Key& key, std::vector<kmrnext::DataPack>& dps,
		   kmrnext::DataStore::MapEnvironment& env)
    {
      {
	assert(dps.size() == 1);
	int nprocs;
	MPI_Comm_size(env.mpi_comm, &nprocs);
	assert(nprocs == 1);
      }
      kmrnext::DataPack& dp = dps.at(0);
      outds->add(key, *dp.data());
      // int val = *(int*)(dp.data()->value());
      // kmrnext::Data d(&val, sizeof(int));
      // outds->add(key, d);
      return 0;
    }
  };

  // A mapper class that collects data in the input DataStores of each
  // process to the master process (keyed-rank) and adds them on the process.
  // Only the master process emits data to the output DataStore.
  class Summarizer0 : public kmrnext::DataStore::Mapper {
    int nprocs_parent_;
    int *masters_;

  public:
    Summarizer0(int nprocs, int *masters) :
      nprocs_parent_(nprocs), masters_(masters) {}
    int operator()(kmrnext::DataStore *inds, kmrnext::DataStore *outds,
		   kmrnext::Key& key, std::vector<kmrnext::DataPack>& dps,
		   kmrnext::DataStore::MapEnvironment& env)
    {
      {
	int nprocs;
	MPI_Comm_size(env.mpi_comm, &nprocs);
	if (nprocs_parent_ >= 4) {
	  assert(dps.size() == 1);
	  assert(nprocs == 4);
	} else {
	  assert(nprocs == nprocs_parent_);
	}
      }
      int master = masters_[key.dim(0)];
      int val = 0;
      if (dps.size() == 1) {
	val = *(int*)(dps.at(0).data()->value());
      } else {
	for (std::vector<kmrnext::DataPack>:: iterator itr = dps.begin();
	     itr != dps.end(); itr++) {
	  val += *(int*)(*itr).data()->value();
	}
      }
      int res;
      MPI_Reduce(&val, &res, 1, MPI_INT, MPI_SUM, master, env.mpi_comm);
      if (env.rank == master) {
	kmrnext::Data d(&res, sizeof(int));
	outds->add(key, d);
      }
      return 0;
    }
  };

  TEST_F(KMRDataStoreTest, Map) {
    // Check owners in case of Serial Mapper
    kmrnext::DataStore ods0(2, gNext);
    ods0.set(ds2_array_);
    Copier0 mapper0;
    kmrnext::View v0(2);
    bool flags0[2] = {true, true};
    v0.set(flags0);
    ds2_->map(&ods0, mapper0, v0);
    EXPECT_EQ(1, *(int*)ods0.get(*k2_00_).data()->value());
    EXPECT_EQ(2, *(int*)ods0.get(*k2_11_).data()->value());
    EXPECT_EQ(3, *(int*)ods0.get(*k2_22_).data()->value());
    EXPECT_EQ(4, *(int*)ods0.get(*k2_33_).data()->value());
    EXPECT_EQ(ds2_owners_[0], ods0.get(*k2_00_).data()->owner());
    EXPECT_EQ(ds2_owners_[1], ods0.get(*k2_11_).data()->owner());
    EXPECT_EQ(ds2_owners_[2], ods0.get(*k2_22_).data()->owner());
    EXPECT_EQ(ds2_owners_[3], ods0.get(*k2_33_).data()->owner());

    // Check owners in case of Parallel Mapper
    kmrnext::DataStore ods1(1, gNext);
    size_t ods1_array[1] = {4};
    ods1.set(ods1_array);
    Summarizer0 mapper1(nprocs, ds2_owners_);
    kmrnext::View v1(2);
    bool flags1[2] = {false, true};
    v1.set(flags1);
    ds2_->map(&ods1, mapper1, v1);
    kmrnext::Key k1_0(1), k1_1(1), k1_2(1), k1_3(1);
    size_t ary_k1_0[1] = {0};
    size_t ary_k1_1[1] = {1};
    size_t ary_k1_2[1] = {2};
    size_t ary_k1_3[1] = {3};
    k1_0.set(ary_k1_0);
    k1_1.set(ary_k1_1);
    k1_2.set(ary_k1_2);
    k1_3.set(ary_k1_3);
    EXPECT_EQ(10, *(int*)ods1.get(k1_0).data()->value());
    EXPECT_EQ(10, *(int*)ods1.get(k1_1).data()->value());
    EXPECT_EQ(10, *(int*)ods1.get(k1_2).data()->value());
    EXPECT_EQ(10, *(int*)ods1.get(k1_3).data()->value());
    EXPECT_EQ(ds2_owners_[0], ods1.get(k1_0).data()->owner());
    EXPECT_EQ(ds2_owners_[1], ods1.get(k1_1).data()->owner());
    EXPECT_EQ(ds2_owners_[2], ods1.get(k1_2).data()->owner());
    EXPECT_EQ(ds2_owners_[3], ods1.get(k1_3).data()->owner());
  }

}
