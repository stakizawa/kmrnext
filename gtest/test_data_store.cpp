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
      key1_ = new Next::Key(ds_size_);
      key2_ = new Next::Key(ds_size_);
      ekey0_ = new Next::Key(ds_size_);
      ekey1_ = new Next::Key(ds_size_);
      size_t ary_key0[3] = {0,0,0};
      size_t ary_key1[3] = {0,1,0};
      size_t ary_key2[3] = {1,1,1};
      size_t ary_ekey0[3] = {0,3,0};
      size_t ary_ekey1[3] = {10,10,10};
      key0_->set(ary_key0);
      key1_->set(ary_key1);
      key2_->set(ary_key2);
      ekey0_->set(ary_ekey0);
      ekey1_->set(ary_ekey1);

      v0_ = new Next::View(ds_size_);
      v1_ = new Next::View(ds_size_);
      v2_ = new Next::View(ds_size_);
      v3_ = new Next::View(ds_size_);
      bool ary_v0[3] = {true, true, true};
      bool ary_v1[3] = {true, false, true};
      bool ary_v2[3] = {true, false, false};
      bool ary_v3[3] = {false, false, false};
      v0_->set(ary_v0);
      v1_->set(ary_v1);
      v2_->set(ary_v2);
      v3_->set(ary_v3);

      d0value_ = 1;
      d0_ = new Next::Data(&d0value_, sizeof(long));
      d1value_ = 2;
      d1_ = new Next::Data(&d1value_, sizeof(long));

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

      ds1_ = new Next::DataStore(2);
      size_t ary_ds1[2] = {2, 2};
      ds1_->set(ary_ds1);
      Next::Key ds1key(2);
      for (size_t i = 0; i < ary_ds1[0]; i++) {
      	ds1key.set_dim(0, i);
      	for (size_t j = 0; j < ary_ds1[1]; j++) {
      	  ds1key.set_dim(1, j);
	  ds1_->add(ds1key, *d1_);
      	}
      }

      ds2_ = new Next::DataStore(2);
      size_t ary_ds2[2] = {2, 1};
      ds2_->set(ary_ds2);
      Next::Key ds2key(2);
      for (size_t i = 0; i < ary_ds2[0]; i++) {
      	ds2key.set_dim(0, i);
      	for (size_t j = 0; j < ary_ds2[1]; j++) {
      	  ds2key.set_dim(1, j);
	  ds2_->add(ds2key, *d1_);
      	}
      }

      ds3_ = new Next::DataStore(1);
      size_t ary_ds3[1] = {2};
      ds3_->set(ary_ds3);
      Next::Key ds3key(1);
      for (size_t i = 0; i < ary_ds3[0]; i++) {
      	ds3key.set_dim(0, i);
	ds3_->add(ds3key, *d1_);
      }
    }

    virtual ~DataStoreTest() {
      delete array_ds0_;
      delete key0_;
      delete key1_;
      delete key2_;
      delete ekey0_;
      delete ekey1_;
      delete v0_;
      delete v1_;
      delete v2_;
      delete v3_;
      delete d0_;
      delete d1_;
      delete ds0_;
      delete ds1_;
      delete ds2_;
      delete ds3_;
    }

    size_t ds_size_;       // 3
    size_t ds_size_error_; // 1000
    size_t *array_ds0_;    // {2,2,2}

    Next::Key *key0_;      // <0,0,0>
    Next::Key *key1_;      // <0,1,0>
    Next::Key *key2_;      // <1,1,1>
    Next::Key *ekey0_;     // <0,3,0>
    Next::Key *ekey1_;     // <10,10,10>

    Next::View *v0_;       // <t,t,t>
    Next::View *v1_;       // <t,f,t>
    Next::View *v2_;       // <t,f,f>
    Next::View *v3_;       // <f,f,f>

    long d0value_;         // 1
    Next::Data *d0_;       // Data(d0value_, sizeof(long))
    long d1value_;         // 2
    Next::Data *d1_;       // Data(d1value_, sizeof(long))

    Next::DataStore *ds0_; // DataStore(Dim:{2,2,2}, Data:d0_)
    Next::DataStore *ds1_; // DataStore(Dim:{2,2},   Data:d1_)
    Next::DataStore *ds2_; // DataStore(Dim:{2,1},   Data:d1_)
    Next::DataStore *ds3_; // DataStore(Dim:{2},     Data:d1_)
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

  TEST_F(DataStoreTest, Get_view) {
    std::vector<Next::DataPack> *vec0 = ds0_->get(*v0_, *key0_);
    EXPECT_EQ((size_t)1, vec0->size());
    delete vec0;
    std::vector<Next::DataPack> *vec1 = ds0_->get(*v1_, *key0_);
    EXPECT_EQ((size_t)2, vec1->size());
    delete vec1;
    std::vector<Next::DataPack> *vec2 = ds0_->get(*v2_, *key0_);
    EXPECT_EQ((size_t)4, vec2->size());
    delete vec2;
    std::vector<Next::DataPack> *vec3 = ds0_->get(*v3_, *key0_);
    EXPECT_EQ((size_t)8, vec3->size());
    delete vec3;

    // If a dimension is out of range, it throws a runtime_error.
    EXPECT_THROW({ds0_->get(*v0_, *ekey0_);}, std::runtime_error);
    // If all dimensions are out of range, it throws a runtime_error.
    EXPECT_THROW({ds0_->get(*v0_, *ekey1_);}, std::runtime_error);
  }

  TEST_F(DataStoreTest, Set_from) {
    std::vector<Next::DataStore*> vec0;
    vec0.push_back(ds1_);
    vec0.push_back(ds1_);
    vec0.push_back(ds1_);
    Next::DataStore mds0(3);
    // assume that ds.get() works fine
    mds0.set_from(vec0);
    EXPECT_EQ((size_t)3, mds0.dim(0));
    EXPECT_EQ(ds1_->dim(0), mds0.dim(1));
    EXPECT_EQ(ds1_->dim(1), mds0.dim(2));
    EXPECT_EQ(*(long*)d1_->value(), *(long*)mds0.get(*key0_).data->value());
    EXPECT_EQ(*(long*)d1_->value(), *(long*)mds0.get(*key1_).data->value());
    EXPECT_EQ(*(long*)d1_->value(), *(long*)mds0.get(*key2_).data->value());

    // If the given vector is empty, it throws a runtime_error.
    std::vector<Next::DataStore*> vec1;
    Next::DataStore mds1(3);
    EXPECT_THROW({mds1.set_from(vec1);}, std::runtime_error);

    // If the given vector has one DataStore, it performs as usual.
    std::vector<Next::DataStore*> vec2;
    vec2.push_back(ds1_);
    Next::DataStore mds2(3);
    mds2.set_from(vec2);
    EXPECT_EQ((size_t)1, mds2.dim(0));
    EXPECT_EQ(ds1_->dim(0), mds2.dim(1));
    EXPECT_EQ(ds1_->dim(1), mds2.dim(2));
    EXPECT_EQ(*(long*)d1_->value(), *(long*)mds2.get(*key0_).data->value());
    EXPECT_EQ(*(long*)d1_->value(), *(long*)mds2.get(*key1_).data->value());
    EXPECT_THROW({mds2.get(*key2_);}, std::runtime_error); // not exist

    // If the data is already set, it throws a runtime_error
    Next::DataStore mds3(3);
    size_t array_mds3[3] = {3, 2, 2};
    mds3.set(array_mds3);
    EXPECT_THROW({mds3.set_from(vec0);}, std::runtime_error);

    // If the dimension sizes of merged DataStores are not same,
    // it throws a runtime_error.
    std::vector<Next::DataStore*> vec4;
    vec4.push_back(ds1_);
    vec4.push_back(ds2_);
    Next::DataStore mds4(3);
    EXPECT_THROW({mds4.set_from(vec4);}, std::runtime_error);

    // If the dimensions of merged DataStores are not same,
    // it throws a runtime_error.
    std::vector<Next::DataStore*> vec5;
    vec5.push_back(ds1_);
    vec5.push_back(ds3_);
    Next::DataStore mds5(3);
    EXPECT_THROW({mds5.set_from(vec5);}, std::runtime_error);
  }

}
