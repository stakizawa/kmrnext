// This test tests DataStore class.
#include <gtest/gtest.h>
#include "kmrnext.hpp"

namespace {

  class DataStoreTest : public ::testing::Test {
  protected:
    DataStoreTest() {
      ds_size_ = 3;
      ds_size_error_ = 1000;
      array_ds0_ = new size_t[ds_size_];
      for (size_t i = 0; i < ds_size_; i++) {
	array_ds0_[i] = 2;
      }

      key0_ = new Next::Key(ds_size_);
      ekey0_ = new Next::Key(ds_size_);
      ekey1_ = new Next::Key(ds_size_);
      size_t ary_key0[3] = {0,0,0};
      size_t ary_ekey0[3] = {0,3,0};
      size_t ary_ekey1[3] = {10,10,10};
      key0_->set(ary_key0);
      ekey0_->set(ary_ekey0);
      ekey1_->set(ary_ekey1);

      d0value_ = 1;
      d0_ = new Next::Data(&d0value_, sizeof(long));

      ds0_ = new Next::DataStore(ds_size_);
      ds0_->set(array_ds0_);
      Next::Key tk(ds_size_);
      for (size_t i = 0; i < 2; i++) {
      	tk.set_dim(0, i);
      	for (size_t j = 0; j < 2; j++) {
      	  tk.set_dim(1, j);
      	  for (size_t k = 0; k < 2; k++) {
      	    tk.set_dim(2, k);
      	    ds0_->add(tk, *d0_);
      	  }
      	}
      }
    }

    virtual ~DataStoreTest() {
      delete array_ds0_;
      delete key0_;
      delete ekey0_;
      delete ekey1_;
      delete d0_;
      delete ds0_;
    }

    virtual void SetUp() {
      // Code here will be called immediately after the constructor (right
      // before each test).
    }

    virtual void TearDown() {
      // Code here will be called immediately after each test (right
      // before the destructor).
    }

    size_t ds_size_;       // 3
    size_t ds_size_error_; // 1000
    size_t *array_ds0_;    // {2,2,2}

    Next::Key *key0_;      // <0,0,0>
    Next::Key *ekey0_;     // <0,3,0>
    Next::Key *ekey1_;     // <10,10,10>

    long d0value_;         // 1
    Next::Data *d0_;       // Data(d0value_, sizeof(long))

    Next::DataStore *ds0_; // DataStore(Dim:{2,2,2}, Data:d0_)
  };

  TEST_F(DataStoreTest, Constructor) {
    Next::DataStore ds(ds_size_);
    EXPECT_EQ(ds_size_, ds.size());
    EXPECT_THROW({Next::DataStore dse(ds_size_error_);}, std::runtime_error);
  }

  TEST_F(DataStoreTest, Set) {
    Next::DataStore ds(ds_size_);
    ds.set(array_ds0_);
    EXPECT_EQ(array_ds0_[0], ds.dim(0));
    EXPECT_EQ(array_ds0_[1], ds.dim(1));
    EXPECT_EQ(array_ds0_[2], ds.dim(2));
    EXPECT_THROW({ds.dim(3);}, std::runtime_error);
    EXPECT_THROW({ds.dim(ds_size_error_);}, std::runtime_error);

    // set can be called only once
    EXPECT_THROW({ds.set(array_ds0_);}, std::runtime_error);
  }

  TEST_F(DataStoreTest, Add) {
    Next::DataStore ds(ds_size_);
    ds.set(array_ds0_);
    // assume that ds.get() works fine
    ds.add(*key0_, *d0_);
    Next::DataPack dp = ds.get(*key0_);
    EXPECT_EQ(*key0_, dp.key);
    EXPECT_EQ(*(long*)d0_->value(), *(long*)dp.data->value());
    EXPECT_EQ(d0_->size(), dp.data->size());

    // If a dimension is out of range, it throws a runtime_error.
    EXPECT_THROW({ds.add(*ekey0_, *d0_);}, std::runtime_error);
    // If all dimensions are out of range, it throws a runtime_error.
    EXPECT_THROW({ds.add(*ekey1_, *d0_);}, std::runtime_error);
  }

  TEST_F(DataStoreTest, Get) {
    Next::DataPack dp = ds0_->get(*key0_);
    EXPECT_EQ(*key0_, dp.key);
    EXPECT_EQ(*(long*)d0_->value(), *(long*)dp.data->value());
    EXPECT_EQ(d0_->size(), dp.data->size());

    // If a dimension is out of range, it throws a runtime_error.
    EXPECT_THROW({ds0_->get(*ekey0_);}, std::runtime_error);
    // If all dimensions are out of range, it throws a runtime_error.
    EXPECT_THROW({ds0_->get(*ekey1_);}, std::runtime_error);
  }
}
