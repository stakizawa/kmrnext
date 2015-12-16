// This test tests DataStore class for KMR backend.
#include <gtest/gtest.h>
#include <sstream>
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

    // void init_data_store(kmrnext::DataStore **dsp, size_t size,
    // 			 size_t *array, kmrnext::Data *d) {
    //   *dsp = new kmrnext::DataStore(size, gNext);
    //   kmrnext::DataStore *ds = *dsp;
    //   ds->set(array);
    //   size_t ndata = 1;
    //   for (size_t i = 0; i < size; i++) {
    // 	ndata *= array[i];
    //   }
    //   for (size_t i = 0; i < ndata; i++) {
    // 	kmrnext::Key k = ds->index_to_key(i);
    // 	ds->add(k, *d);
    //   }
    // }

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

#if 0
  TEST_F(KMRDataStoreTest, Get_view) {
    // dummy codes
    EXPECT_EQ(1, *(int*)ds2_->get(*k2_00_).data()->value());
    EXPECT_EQ(2, *(int*)ds2_->get(*k2_11_).data()->value());
    EXPECT_EQ(3, *(int*)ds2_->get(*k2_22_).data()->value());
    EXPECT_EQ(4, *(int*)ds2_->get(*k2_33_).data()->value());
    EXPECT_EQ(ds2_owners_[0], ds2_->get(*k2_00_).data()->owner());
    EXPECT_EQ(ds2_owners_[1], ds2_->get(*k2_11_).data()->owner());
    EXPECT_EQ(ds2_owners_[2], ds2_->get(*k2_22_).data()->owner());
    EXPECT_EQ(ds2_owners_[3], ds2_->get(*k2_33_).data()->owner());
  }
#endif

#if 0
  TEST_F(KMRDataStoreTest, Set_from) {
  }
#endif

#if 0
  TEST_F(KMRDataStoreTest, Split_to) {
  }
#endif

#if 0
  TEST_F(KMRDataStoreTest, Map) {
  }
#endif

}
