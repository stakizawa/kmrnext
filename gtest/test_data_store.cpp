// This test tests DataStore class.
#include <gtest/gtest.h>
#include <sstream>
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

      key2d0_ = new Next::Key(2);
      key2d1_ = new Next::Key(2);
      size_t ary_key2d0[2] = {0,0};
      size_t ary_key2d1[2] = {1,0};
      key2d0_->set(ary_key2d0);
      key2d1_->set(ary_key2d1);

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
      delete key2d0_;
      delete key2d1_;
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

    Next::Key *key2d0_;    // <0,0>
    Next::Key *key2d1_;    // <1,0>

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
    // assume that ds.get() works fine
    Next::DataStore ds(ds_size_);
    ds.set(array_ds0_);
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

    // If dimension sizes of Key and DataStore are not same,
    // it throws a runtime_error.
    EXPECT_THROW({ds0_->get(*key2d0_);}, std::runtime_error);
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
    // assume that ds.get() works fine
    std::vector<Next::DataStore*> vec0;
    vec0.push_back(ds1_);
    vec0.push_back(ds1_);
    vec0.push_back(ds1_);
    Next::DataStore mds0(3);
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

  TEST_F(DataStoreTest, Split_to) {
    // assume that ds.get() works fine
    std::vector<Next::DataStore*> vec0;
    Next::DataStore sds00(ds_size_ - 1);
    Next::DataStore sds01(ds_size_ - 1);
    vec0.push_back(&sds00);
    vec0.push_back(&sds01);
    ds0_->split_to(vec0);
    EXPECT_EQ((size_t)2, vec0.size());
    EXPECT_EQ(array_ds0_[1], sds00.dim(0));
    EXPECT_EQ(array_ds0_[2], sds00.dim(1));
    EXPECT_EQ(*(long*)d0_->value(), *(long*)sds00.get(*key2d0_).data->value());
    EXPECT_EQ(*(long*)d0_->value(), *(long*)sds00.get(*key2d1_).data->value());
    EXPECT_EQ(array_ds0_[1], sds01.dim(0));
    EXPECT_EQ(array_ds0_[2], sds01.dim(1));
    EXPECT_EQ(*(long*)d0_->value(), *(long*)sds01.get(*key2d0_).data->value());
    EXPECT_EQ(*(long*)d0_->value(), *(long*)sds01.get(*key2d1_).data->value());

    // If target DataStore has only one dimension, it throws runtime_errror.
    std::vector<Next::DataStore*> vec1;
    Next::DataStore sds10(1);
    Next::DataStore sds11(1);
    vec1.push_back(&sds10);
    vec1.push_back(&sds11);
    EXPECT_THROW({ds3_->split_to(vec1);}, std::runtime_error);

    // If the vector is empty, it throws runtime_error.
    std::vector<Next::DataStore*> vec2;
    EXPECT_THROW({ds0_->split_to(vec2);}, std::runtime_error);

    // If more vector than the top-level dimension of the target DataStore
    // is given, it throws runtime_error.
    std::vector<Next::DataStore*> vec3;
    Next::DataStore sds30(ds_size_ - 1);
    Next::DataStore sds31(ds_size_ - 1);
    Next::DataStore sds32(ds_size_ - 1);
    vec3.push_back(&sds30);
    vec3.push_back(&sds31);
    vec3.push_back(&sds32);
    EXPECT_THROW({ds0_->split_to(vec3);}, std::runtime_error);
  }

  // A mapper class that increments value and the calculates average.
  class Summarizer {
  public:
    int operator()(Next::DataStore *inds, Next::DataStore *outds,
		   Next::Key key, std::vector<Next::DataPack>& dps)
    {
      long sum = 0;
      for (size_t i = 0; i < dps.size(); i++) {
	Next::DataPack& dp = dps.at(i);
	long v = *(long *)dp.data->value() + 1;
	sum += v;
      }
      long avg = sum / dps.size();
      Next::Data d(&avg, sizeof(long));
      outds->add(key, d);
      return 0;
    }
  };

  TEST_F(DataStoreTest, Map) {
    // assume that ds.get() works fine
    Summarizer mapper;

    Next::DataStore ods0(ds_size_);
    ods0.set(array_ds0_);
    ds0_->map(&ods0, mapper, *v0_);
    EXPECT_EQ(2, *(long*)ods0.get(*key0_).data->value());
    EXPECT_EQ(2, *(long*)ods0.get(*key1_).data->value());
    EXPECT_EQ(2, *(long*)ods0.get(*key2_).data->value());

    Next::DataStore ods1(2);
    size_t ary_ods1[2] = {2,2};
    ods1.set(ary_ods1);
    ds0_->map(&ods1, mapper, *v1_);
    EXPECT_EQ(2, *(long*)ods1.get(*key2d0_).data->value());
    EXPECT_EQ(2, *(long*)ods1.get(*key2d1_).data->value());

    // If the input and output DataStores are same, it throws runtime_error.
    EXPECT_THROW({ds0_->map(ds0_, mapper, *v0_);}, std::runtime_error);

    // If the output DataStore already has some value, it throws runtime_error.
    // TODO check future.

    // If the dimension of view does not match that of the input DataStore,
    // it throws runtime_error.
    Next::DataStore ods2(ds_size_);
    ods2.set(array_ds0_);
    EXPECT_THROW({ds1_->map(&ods2, mapper, *v0_);}, std::runtime_error);
  }

  class DataLoader1D {
    size_t size_;
  public:
    DataLoader1D(size_t siz) : size_(siz) {}

    int operator()(Next::DataStore *ds, int num)
    {
      Next::Key key(1);
      Next::Data data(&num, sizeof(int));
      for (size_t i = 0; i < size_; i++) {
	key.set_dim(0, i);
	ds->add(key, data);
      }
      return 0;
    }
  };

  class DataLoader2D {
    size_t size0_;
    size_t size1_;
  public:
    DataLoader2D(size_t siz0, size_t siz1) : size0_(siz0), size1_(siz1) {}

    int operator()(Next::DataStore *ds, int num)
    {
      Next::Key key(2);
      Next::Data data(&num, sizeof(int));
      for (size_t i = 0; i < size0_; i++) {
	key.set_dim(0, i);
	for (size_t j = 0; j < size1_; j++) {
	  key.set_dim(1, j);
	  ds->add(key, data);
	}
      }
      return 0;
    }
  };

  TEST_F(DataStoreTest, Load_array) {
    DataLoader1D loader1d(ds0_->dim(2));
    DataLoader2D loader2d(ds0_->dim(1), ds0_->dim(2));

    // Test 0
    Next::DataStore ds0(ds_size_);
    ds0.set(array_ds0_);
    std::vector<int> vec0;
    for (size_t i = 0; i < ds0.dim(0); i++) {
      for (size_t j = 0; j < ds0.dim(1); j++) {
	vec0.push_back(j+1);
      }
    }
    ds0.load_array(vec0, loader1d);
    EXPECT_EQ(1, *(long*)ds0.get(*key0_).data->value());
    EXPECT_EQ(2, *(long*)ds0.get(*key1_).data->value());
    EXPECT_EQ(2, *(long*)ds0.get(*key2_).data->value());

    // Test 1
    Next::DataStore ds1(ds_size_);
    ds1.set(array_ds0_);
    std::vector<int> vec1;
    for (size_t i = 0; i < ds1.dim(0); i++) {
      vec1.push_back(i+1);
    }
    ds1.load_array(vec1, loader2d);
    EXPECT_EQ(1, *(long*)ds1.get(*key0_).data->value());
    EXPECT_EQ(1, *(long*)ds1.get(*key1_).data->value());
    EXPECT_EQ(2, *(long*)ds1.get(*key2_).data->value());

    // IF the size of array is not same as the multiple of dimension size
    // of the DataStore, it throws runtime_error.
    Next::DataStore ds2(ds_size_);
    ds2.set(array_ds0_);
    std::vector<int> vec2;
    for (size_t i = 0; i < ds_size_error_; i++) {
      vec2.push_back(i+1);
    }
    EXPECT_THROW({ds2.load_array(vec2, loader2d);}, std::runtime_error);
  }

  class DS0Printer : public Next::DataPackDumper {
  public:
    std::string operator()(Next::DataPack& dp)
    {
      std::ostringstream os;
      os << *(long*)dp.data->value() << ",";
      return os.str();
    }
  };

  TEST_F(DataStoreTest, Dump) {
    DS0Printer ptr0;
    std::string expected0 = "Data Count: 8\n1,1,1,1,1,1,1,1,";
    std::string actual0 = ds0_->dump(ptr0);
    EXPECT_STREQ(expected0.c_str(), actual0.c_str());
  }

}
