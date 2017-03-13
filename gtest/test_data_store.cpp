// This test tests DataStore class.
#include <gtest/gtest.h>
#include <sstream>
#include "kmrnext.hpp"

extern kmrnext::KMRNext *gNext;

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

      key0_ = new kmrnext::Key(ds_size_);
      key1_ = new kmrnext::Key(ds_size_);
      key2_ = new kmrnext::Key(ds_size_);
      ekey0_ = new kmrnext::Key(ds_size_);
      ekey1_ = new kmrnext::Key(ds_size_);
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

      key2d0_ = new kmrnext::Key(2);
      key2d1_ = new kmrnext::Key(2);
      size_t ary_key2d0[2] = {0,0};
      size_t ary_key2d1[2] = {1,0};
      key2d0_->set(ary_key2d0);
      key2d1_->set(ary_key2d1);

      v0_ = new kmrnext::View(ds_size_);
      v1_ = new kmrnext::View(ds_size_);
      v2_ = new kmrnext::View(ds_size_);
      v3_ = new kmrnext::View(ds_size_);
      long ary_v0[3] = { kmrnext::View::SplitAll,
			 kmrnext::View::SplitAll,
			 kmrnext::View::SplitAll };
      long ary_v1[3] = { kmrnext::View::SplitAll,
			 kmrnext::View::SplitNone,
			 kmrnext::View::SplitAll };
      long ary_v2[3] = { kmrnext::View::SplitAll,
			 kmrnext::View::SplitNone,
			 kmrnext::View::SplitNone };
      long ary_v3[3] = { kmrnext::View::SplitNone,
			 kmrnext::View::SplitNone,
			 kmrnext::View::SplitNone };
      v0_->set(ary_v0);
      v1_->set(ary_v1);
      v2_->set(ary_v2);
      v3_->set(ary_v3);

      d0value_ = 1;
      d0_ = new kmrnext::Data(&d0value_, sizeof(long));
      d1value_ = 2;
      d1_ = new kmrnext::Data(&d1value_, sizeof(long));

      ds0_ = gNext->create_ds(ds_size_);
      ds0_->set(array_ds0_);
      kmrnext::Key tk(ds_size_);
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

      ds1_ = gNext->create_ds(2);
      size_t ary_ds1[2] = {2, 2};
      ds1_->set(ary_ds1);
      kmrnext::Key ds1key(2);
      for (size_t i = 0; i < ary_ds1[0]; i++) {
      	ds1key.set_dim(0, i);
      	for (size_t j = 0; j < ary_ds1[1]; j++) {
      	  ds1key.set_dim(1, j);
	  ds1_->add(ds1key, *d1_);
      	}
      }

      ds2_ = gNext->create_ds(2);
      size_t ary_ds2[2] = {2, 1};
      ds2_->set(ary_ds2);
      kmrnext::Key ds2key(2);
      for (size_t i = 0; i < ary_ds2[0]; i++) {
      	ds2key.set_dim(0, i);
      	for (size_t j = 0; j < ary_ds2[1]; j++) {
      	  ds2key.set_dim(1, j);
	  ds2_->add(ds2key, *d1_);
      	}
      }

      ds3_ = gNext->create_ds(1);
      size_t ary_ds3[1] = {2};
      ds3_->set(ary_ds3);
      kmrnext::Key ds3key(1);
      for (size_t i = 0; i < ary_ds3[0]; i++) {
      	ds3key.set_dim(0, i);
	ds3_->add(ds3key, *d1_);
      }

      ds4_ = gNext->create_ds(3);
      ds4_->set(array_ds0_);
    }

    virtual ~DataStoreTest() {
      delete[] array_ds0_;
      delete   key0_;
      delete   key1_;
      delete   key2_;
      delete   ekey0_;
      delete   ekey1_;
      delete   key2d0_;
      delete   key2d1_;
      delete   v0_;
      delete   v1_;
      delete   v2_;
      delete   v3_;
      delete   d0_;
      delete   d1_;
      delete   ds0_;
      delete   ds1_;
      delete   ds2_;
      delete   ds3_;
      delete   ds4_;
    }

    size_t ds_size_;          // 3
    size_t ds_size_error_;    // 1000
    size_t *array_ds0_;       // {2,2,2}

    kmrnext::Key *key0_;      // <0,0,0>
    kmrnext::Key *key1_;      // <0,1,0>
    kmrnext::Key *key2_;      // <1,1,1>
    kmrnext::Key *ekey0_;     // <0,3,0>
    kmrnext::Key *ekey1_;     // <10,10,10>

    kmrnext::Key *key2d0_;    // <0,0>
    kmrnext::Key *key2d1_;    // <1,0>

    kmrnext::View *v0_;       // <t,t,t>
    kmrnext::View *v1_;       // <t,f,t>
    kmrnext::View *v2_;       // <t,f,f>
    kmrnext::View *v3_;       // <f,f,f>

    long d0value_;            // 1
    kmrnext::Data *d0_;       // Data(d0value_, sizeof(long))
    long d1value_;            // 2
    kmrnext::Data *d1_;       // Data(d1value_, sizeof(long))

    kmrnext::DataStore *ds0_; // DataStore(Dim:{2,2,2}, Data:d0_)
    kmrnext::DataStore *ds1_; // DataStore(Dim:{2,2},   Data:d1_)
    kmrnext::DataStore *ds2_; // DataStore(Dim:{2,1},   Data:d1_)
    kmrnext::DataStore *ds3_; // DataStore(Dim:{2},     Data:d1_)
    kmrnext::DataStore *ds4_; // DataStore(Dim:{2,2,2}, Data:NULL)
  };

  TEST_F(DataStoreTest, Constructor) {
    kmrnext::DataStore* ds = gNext->create_ds(ds_size_);
    EXPECT_EQ(ds_size_, ds->size());
    EXPECT_THROW({ gNext->create_ds(ds_size_error_); },
		 std::runtime_error);
  }

  TEST_F(DataStoreTest, Set) {
    kmrnext::DataStore *ds = gNext->create_ds(ds_size_);
    ds->set(array_ds0_);
    EXPECT_EQ(array_ds0_[0], ds->dim(0));
    EXPECT_EQ(array_ds0_[1], ds->dim(1));
    EXPECT_EQ(array_ds0_[2], ds->dim(2));
    EXPECT_THROW({ds->dim(3);}, std::runtime_error);
    EXPECT_THROW({ds->dim(ds_size_error_);}, std::runtime_error);

    // set can be called only once
    EXPECT_THROW({ds->set(array_ds0_);}, std::runtime_error);
    delete ds;
  }

  TEST_F(DataStoreTest, Add) {
    // assume that ds.get() works fine
    kmrnext::DataStore *ds = gNext->create_ds(ds_size_);
    ds->set(array_ds0_);
    ds->add(*key0_, *d0_);
    kmrnext::DataPack dp = ds->get(*key0_);
    EXPECT_EQ(*key0_, dp.key());
    EXPECT_EQ(*static_cast<long*>(d0_->value()),
	      *static_cast<long*>(dp.data().value()));
    EXPECT_EQ(d0_->size(), dp.data().size());

    // If a dimension is out of range, it throws a runtime_error.
    EXPECT_THROW({ds->add(*ekey0_, *d0_);}, std::runtime_error);
    // If all dimensions are out of range, it throws a runtime_error.
    EXPECT_THROW({ds->add(*ekey1_, *d0_);}, std::runtime_error);
    delete ds;
  }

  TEST_F(DataStoreTest, Get) {
    kmrnext::DataPack dp0 = ds0_->get(*key0_);
    EXPECT_EQ(*key0_, dp0.key());
    EXPECT_EQ(*static_cast<long*>(d0_->value()),
	      *static_cast<long*>(dp0.data().value()));
    EXPECT_EQ(d0_->size(), dp0.data().size());

    // If the data of the specified Key is not set, it returns NULL.
    kmrnext::DataPack dp1 = ds4_->get(*key0_);
    EXPECT_EQ(*key0_, dp1.key());
    EXPECT_EQ(NULL, dp1.data().value());

    // If dimension sizes of Key and DataStore are not same,
    // it throws a runtime_error.
    EXPECT_THROW({ds0_->get(*key2d0_);}, std::runtime_error);
    // If a dimension is out of range, it throws a runtime_error.
    EXPECT_THROW({ds0_->get(*ekey0_);}, std::runtime_error);
    // If all dimensions are out of range, it throws a runtime_error.
    EXPECT_THROW({ds0_->get(*ekey1_);}, std::runtime_error);
  }

  TEST_F(DataStoreTest, Get_view) {
    std::vector<kmrnext::DataPack> *vec0 = ds0_->get(*v0_, *key0_);
    EXPECT_EQ(1, vec0->size());
    delete vec0;
    std::vector<kmrnext::DataPack> *vec1 = ds0_->get(*v1_, *key0_);
    EXPECT_EQ(2, vec1->size());
    delete vec1;
    std::vector<kmrnext::DataPack> *vec2 = ds0_->get(*v2_, *key0_);
    EXPECT_EQ(4, vec2->size());
    delete vec2;
    std::vector<kmrnext::DataPack> *vec3 = ds0_->get(*v3_, *key0_);
    EXPECT_EQ(8, vec3->size());
    delete vec3;

    // If a dimension is out of range, it throws a runtime_error.
    EXPECT_THROW({ds0_->get(*v0_, *ekey0_);}, std::runtime_error);
    // If all dimensions are out of range, it throws a runtime_error.
    EXPECT_THROW({ds0_->get(*v0_, *ekey1_);}, std::runtime_error);
  }

  TEST_F(DataStoreTest, Remove) {
    kmrnext::DataPack dp0 = ds0_->remove(*key0_);
    EXPECT_EQ(*key0_, dp0.key());
    EXPECT_EQ(*static_cast<long*>(d0_->value()),
	      *static_cast<long*>(dp0.data().value()));
    EXPECT_EQ(d0_->size(), dp0.data().size());
    // If the data is correctly remove, the following get returns NULL.
    kmrnext::DataPack dp1 = ds0_->get(*key0_);
    EXPECT_EQ(*key0_, dp1.key());
    EXPECT_EQ(NULL, dp1.data().value());
    // If the data is correctly remove, the new Data can be set to
    // the same Key.
    EXPECT_NO_THROW({ds0_->add(*key0_, *d1_);});
    kmrnext::DataPack dp2 = ds0_->get(*key0_);
    EXPECT_EQ(*key0_, dp2.key());
    EXPECT_EQ(*static_cast<long*>(d1_->value()),
	      *static_cast<long*>(dp2.data().value()));
    EXPECT_EQ(d1_->size(), dp2.data().size());

    // If the data of the specified Kye is not set, it returns NULL.
    kmrnext::DataPack dp3 = ds4_->remove(*key0_);
    EXPECT_EQ(*key0_, dp3.key());
    EXPECT_EQ(NULL, dp3.data().value());

    // If dimension sizes of Key and DataStore are not same,
    // it throws a runtime_error.
    EXPECT_THROW({ds0_->get(*key2d0_);}, std::runtime_error);
    // If a dimension is out of range, it throws a runtime_error.
    EXPECT_THROW({ds0_->get(*ekey0_);}, std::runtime_error);
    // If all dimensions are out of range, it throws a runtime_error.
    EXPECT_THROW({ds0_->get(*ekey1_);}, std::runtime_error);
  }

  TEST_F(DataStoreTest, Set_from) {
    // assume that ds.get() works fine
    std::vector<kmrnext::DataStore*> vec0;
    vec0.push_back(ds1_);
    vec0.push_back(ds1_);
    vec0.push_back(ds1_);
    kmrnext::DataStore *mds0 = gNext->create_ds(3);
    mds0->set_from(vec0);
    EXPECT_EQ(3, mds0->dim(0));
    EXPECT_EQ(ds1_->dim(0), mds0->dim(1));
    EXPECT_EQ(ds1_->dim(1), mds0->dim(2));
    EXPECT_EQ(*static_cast<long*>(d1_->value()),
	      *static_cast<long*>(mds0->get(*key0_).data().value()));
    EXPECT_EQ(*static_cast<long*>(d1_->value()),
	      *static_cast<long*>(mds0->get(*key1_).data().value()));
    EXPECT_EQ(*static_cast<long*>(d1_->value()),
	      *static_cast<long*>(mds0->get(*key2_).data().value()));
    delete mds0;

    // If the given vector is empty, it throws a runtime_error.
    std::vector<kmrnext::DataStore*> vec1;
    kmrnext::DataStore *mds1 = gNext->create_ds(3);
    EXPECT_THROW({mds1->set_from(vec1);}, std::runtime_error);
    delete mds1;

    // If the given vector has one DataStore, it performs as usual.
    std::vector<kmrnext::DataStore*> vec2;
    vec2.push_back(ds1_);
    kmrnext::DataStore *mds2 = gNext->create_ds(3);
    mds2->set_from(vec2);
    EXPECT_EQ(1, mds2->dim(0));
    EXPECT_EQ(ds1_->dim(0), mds2->dim(1));
    EXPECT_EQ(ds1_->dim(1), mds2->dim(2));
    EXPECT_EQ(*static_cast<long*>(d1_->value()),
	      *static_cast<long*>(mds2->get(*key0_).data().value()));
    EXPECT_EQ(*static_cast<long*>(d1_->value()),
	      *static_cast<long*>(mds2->get(*key1_).data().value()));
    EXPECT_THROW({mds2->get(*key2_);}, std::runtime_error); // not exist
    delete mds2;

    // If the data is already set, it throws a runtime_error
    kmrnext::DataStore *mds3 = gNext->create_ds(3);
    size_t array_mds3[3] = {3, 2, 2};
    mds3->set(array_mds3);
    EXPECT_THROW({mds3->set_from(vec0);}, std::runtime_error);
    delete mds3;

    // If the dimension sizes of merged DataStores are not same,
    // it throws a runtime_error.
    std::vector<kmrnext::DataStore*> vec4;
    vec4.push_back(ds1_);
    vec4.push_back(ds2_);
    kmrnext::DataStore *mds4 = gNext->create_ds(3);
    EXPECT_THROW({mds4->set_from(vec4);}, std::runtime_error);
    delete mds4;

    // If the dimensions of merged DataStores are not same,
    // it throws a runtime_error.
    std::vector<kmrnext::DataStore*> vec5;
    vec5.push_back(ds1_);
    vec5.push_back(ds3_);
    kmrnext::DataStore *mds5 = gNext->create_ds(3);
    EXPECT_THROW({mds5->set_from(vec5);}, std::runtime_error);
    delete mds5;
  }

  TEST_F(DataStoreTest, Split_to) {
    // assume that ds.get() works fine
    std::vector<kmrnext::DataStore*> vec0;
    kmrnext::DataStore *sds00 = gNext->create_ds(ds_size_ - 1);
    kmrnext::DataStore *sds01 = gNext->create_ds(ds_size_ - 1);
    vec0.push_back(sds00);
    vec0.push_back(sds01);
    ds0_->split_to(vec0);
    EXPECT_EQ(2, vec0.size());
    EXPECT_EQ(array_ds0_[1], sds00->dim(0));
    EXPECT_EQ(array_ds0_[2], sds00->dim(1));
    EXPECT_EQ(*static_cast<long*>(d0_->value()),
	      *static_cast<long*>(sds00->get(*key2d0_).data().value()));
    EXPECT_EQ(*static_cast<long*>(d0_->value()),
	      *static_cast<long*>(sds00->get(*key2d1_).data().value()));
    EXPECT_EQ(array_ds0_[1], sds01->dim(0));
    EXPECT_EQ(array_ds0_[2], sds01->dim(1));
    EXPECT_EQ(*static_cast<long*>(d0_->value()),
	      *static_cast<long*>(sds01->get(*key2d0_).data().value()));
    EXPECT_EQ(*static_cast<long*>(d0_->value()),
	      *static_cast<long*>(sds01->get(*key2d1_).data().value()));
    delete sds00;
    delete sds01;

    // If target DataStore has only one dimension, it throws runtime_errror.
    std::vector<kmrnext::DataStore*> vec1;
    kmrnext::DataStore *sds10 = gNext->create_ds(1);
    kmrnext::DataStore *sds11 = gNext->create_ds(1);
    vec1.push_back(sds10);
    vec1.push_back(sds11);
    EXPECT_THROW({ds3_->split_to(vec1);}, std::runtime_error);
    delete sds10;
    delete sds11;

    // If the vector is empty, it throws runtime_error.
    std::vector<kmrnext::DataStore*> vec2;
    EXPECT_THROW({ds0_->split_to(vec2);}, std::runtime_error);

    // If more vector than the top-level dimension of the target DataStore
    // is given, it throws runtime_error.
    std::vector<kmrnext::DataStore*> vec3;
    kmrnext::DataStore *sds30 = gNext->create_ds(ds_size_ - 1);
    kmrnext::DataStore *sds31 = gNext->create_ds(ds_size_ - 1);
    kmrnext::DataStore *sds32 = gNext->create_ds(ds_size_ - 1);
    vec3.push_back(sds30);
    vec3.push_back(sds31);
    vec3.push_back(sds32);
    EXPECT_THROW({ds0_->split_to(vec3);}, std::runtime_error);
    delete sds30;
    delete sds31;
    delete sds32;
  }

  // A mapper class that increments value and the calculates average.
  class Summarizer : public kmrnext::DataStore::Mapper {
  public:
    int operator()(kmrnext::DataStore *inds, kmrnext::DataStore *outds,
		   kmrnext::Key& key, std::vector<kmrnext::DataPack>& dps,
		   kmrnext::DataStore::MapEnvironment& env)
    {
      long sum = 0;
      for (size_t i = 0; i < dps.size(); i++) {
	kmrnext::DataPack& dp = dps.at(i);
	long v = *static_cast<long*>(dp.data().value()) + 1;
	sum += v;
      }
      long avg = sum / static_cast<long>(dps.size());
      kmrnext::Data d(&avg, sizeof(long));
      outds->add(key, d);
      return 0;
    }
  };

  // A mapper class that always fails.
  class FailMapper : public kmrnext::DataStore::Mapper {
  public:
    int operator()(kmrnext::DataStore *inds, kmrnext::DataStore *outds,
		   kmrnext::Key& key, std::vector<kmrnext::DataPack>& dps,
		   kmrnext::DataStore::MapEnvironment& env)
    {
      return 1;
    }
  };

  TEST_F(DataStoreTest, Map) {
    // assume that ds.get() works fine
    Summarizer mapper;

    kmrnext::DataStore *ods0 = gNext->create_ds(ds_size_);
    ods0->set(array_ds0_);
    ds0_->map(mapper, *v0_, ods0);
    EXPECT_EQ(2, *static_cast<long*>(ods0->get(*key0_).data().value()));
    EXPECT_EQ(2, *static_cast<long*>(ods0->get(*key1_).data().value()));
    EXPECT_EQ(2, *static_cast<long*>(ods0->get(*key2_).data().value()));
    delete ods0;

    kmrnext::DataStore *ods1 = gNext->create_ds(2);
    size_t ary_ods1[2] = {2,2};
    ods1->set(ary_ods1);
    ds0_->map(mapper, *v1_, ods1);
    EXPECT_EQ(2, *static_cast<long*>(ods1->get(*key2d0_).data().value()));
    EXPECT_EQ(2, *static_cast<long*>(ods1->get(*key2d1_).data().value()));
    delete ods1;

    // If the dimension of view does not match that of the input DataStore,
    // it throws runtime_error.
    kmrnext::DataStore *ods2 = gNext->create_ds(ds_size_);
    ods2->set(array_ds0_);
    EXPECT_THROW({ds1_->map(mapper, *v0_, ods2);}, std::runtime_error);
    delete ods2;

    // If NULL is given to the output DataStore, it throws runtime_error.
    EXPECT_THROW({ds0_->map(mapper, *v0_, NULL);}, std::runtime_error);

    // If the output DataStore is omitted, map overwrites the input DataStore.
    kmrnext::DataStore *ds0 = ds0_->duplicate();
    EXPECT_EQ(1, *static_cast<long*>(ds0->get(*key0_).data().value()));
    ds0->map(mapper, *v0_);
    EXPECT_EQ(2, *static_cast<long*>(ds0->get(*key0_).data().value()));
    delete ds0;

#ifdef BACKEND_SERIAL
    // This test destroys KMRNext environment.
    // If a map task returns a number other than 0, it throws runtime_error.
    FailMapper failer;
    ds0 = ds0_->duplicate();
    EXPECT_THROW({ds0->map(failer, *v0_);}, std::runtime_error);
    delete ds0;
#endif
  }

  class DataLoader1D : public kmrnext::DataStore::Loader<long> {
    size_t key_size_;
    size_t data_count_;
  public:
    DataLoader1D(size_t key_siz, size_t dat_cnt)
      : key_size_(key_siz), data_count_(dat_cnt) {}

    int operator()(kmrnext::DataStore *ds, const long& num)
    {
      kmrnext::Key key(key_size_);
      size_t k0 = static_cast<size_t>(num) / ds->dim(1);
      size_t k1 = static_cast<size_t>(num) % ds->dim(1);
      key.set_dim(0, k0);
      key.set_dim(1, k1);
      long val = static_cast<long>(k1) + 1;
      kmrnext::Data data(static_cast<void*>(&val), sizeof(long));
      for (size_t i = 0; i < data_count_; i++) {
	key.set_dim(2, i);
	ds->add(key, data);
      }
      return 0;
    }
  };

  class DataLoader2D : public kmrnext::DataStore::Loader<long> {
    size_t key_size_;
    size_t data_count0_;
    size_t data_count1_;
  public:
    DataLoader2D(size_t key_siz, size_t cnt0, size_t cnt1)
      : key_size_(key_siz), data_count0_(cnt0), data_count1_(cnt1) {}

    int operator()(kmrnext::DataStore *ds, const long& num)
    {
      kmrnext::Key key(key_size_);
      key.set_dim(0, static_cast<size_t>(num));
      long val = num + 1;
      kmrnext::Data data(static_cast<void*>(&val), sizeof(long));
      for (size_t i = 0; i < data_count0_; i++) {
	key.set_dim(1, i);
	for (size_t j = 0; j < data_count1_; j++) {
	  key.set_dim(2, j);
	  ds->add(key, data);
	}
      }
      return 0;
    }
  };

  class DataLoader3D : public kmrnext::DataStore::Loader<long> {
    size_t key_size_;
    size_t data_count0_;
    size_t data_count1_;
    size_t data_count2_;
  public:
    DataLoader3D(size_t key_siz, size_t cnt0, size_t cnt1, size_t cnt2)
      : key_size_(key_siz), data_count0_(cnt0), data_count1_(cnt1),
	data_count2_(cnt2){}

    int operator()(kmrnext::DataStore *ds, const long& num)
    {
      kmrnext::Key key(key_size_);
      long val = num + 1;
      kmrnext::Data data(static_cast<void*>(&val), sizeof(long));
      for (size_t i = 0; i < data_count0_; i++) {
	key.set_dim(0, i);
	for (size_t j = 0; j < data_count1_; j++) {
	  key.set_dim(1, j);
	  for (size_t k = 0; k < data_count2_; k++) {
	    key.set_dim(2, k);
	    ds->add(key, data);
	  }
	}
      }
      return 0;
    }
  };

  TEST_F(DataStoreTest, Load_array) {
    DataLoader1D loader1d(ds_size_, ds0_->dim(2));
    DataLoader2D loader2d(ds_size_, ds0_->dim(1), ds0_->dim(2));
    DataLoader3D loader3d(ds_size_, ds0_->dim(0), ds0_->dim(1), ds0_->dim(2));

    // Test 0: DataLoader1D, load 2x2 times
    kmrnext::DataStore *ds0 = gNext->create_ds(ds_size_);
    ds0->set(array_ds0_);
    std::vector<long> vec0;
    for (size_t i = 0; i < ds0->dim(0); i++) {
      for (size_t j = 0; j < ds0->dim(1); j++) {
	vec0.push_back(static_cast<long>(i * ds0->dim(1) + j));
      }
    }
    ds0->load_integers(vec0, loader1d);
    EXPECT_EQ(1, *static_cast<long*>(ds0->get(*key0_).data().value()));
    EXPECT_EQ(2, *static_cast<long*>(ds0->get(*key1_).data().value()));
    EXPECT_EQ(2, *static_cast<long*>(ds0->get(*key2_).data().value()));
    delete ds0;

    // Test 1: DataLoader2D, load 2 times
    kmrnext::DataStore *ds1 = gNext->create_ds(ds_size_);
    ds1->set(array_ds0_);
    std::vector<long> vec1;
    for (size_t i = 0; i < ds1->dim(0); i++) {
      vec1.push_back(static_cast<long>(i));
    }
    ds1->load_integers(vec1, loader2d);
    EXPECT_EQ(1, *static_cast<long*>(ds1->get(*key0_).data().value()));
    EXPECT_EQ(1, *static_cast<long*>(ds1->get(*key1_).data().value()));
    EXPECT_EQ(2, *static_cast<long*>(ds1->get(*key2_).data().value()));
    delete ds1;

    // Test 2: DataLoader3D, load once
    kmrnext::DataStore *ds2 = gNext->create_ds(ds_size_);
    ds2->set(array_ds0_);
    std::vector<long> vec2;
    vec2.push_back(0);
    ds2->load_integers(vec2, loader3d);
    EXPECT_EQ(1, *static_cast<long*>(ds2->get(*key0_).data().value()));
    EXPECT_EQ(1, *static_cast<long*>(ds2->get(*key1_).data().value()));
    EXPECT_EQ(1, *static_cast<long*>(ds2->get(*key2_).data().value()));
    delete ds2;

    // If the size of array is not same as the product of dimension sizes
    // of the DataStore, it throws runtime_error.
    kmrnext::DataStore *ds3 = gNext->create_ds(ds_size_);
    ds3->set(array_ds0_);
    std::vector<long> vec3;
    for (size_t i = 0; i < ds_size_error_; i++) {
      vec3.push_back(static_cast<long>(i));
    }
    EXPECT_THROW({ds3->load_integers(vec3, loader2d);}, std::runtime_error);
    delete ds3;
  }

  class DS0Printer : public kmrnext::DataPack::Dumper {
  public:
    std::string operator()(kmrnext::DataPack& dp)
    {
      std::ostringstream os;
      os << *static_cast<long*>(dp.data().value()) << ",";
      return os.str();
    }
  };

  TEST_F(DataStoreTest, Dump) {
    DS0Printer ptr0;
    std::string expected0 = "1,1,1,1,1,1,1,1,";
    std::string actual0 = ds0_->dump(ptr0);
    EXPECT_STREQ(expected0.c_str(), actual0.c_str());

    // If the DataStore is not inisialized yet, it returns "".
    kmrnext::DataStore *ds1 = gNext->create_ds(ds_size_);
    std::string expected1 = "";
    std::string actual1 = ds1->dump(ptr0);
    EXPECT_STREQ(expected1.c_str(), actual1.c_str());
    ds1->set(array_ds0_);
    std::string expected2 = "";
    std::string actual2 = ds1->dump(ptr0);
    EXPECT_STREQ(expected2.c_str(), actual2.c_str());
    delete ds1;
  }

  TEST_F(DataStoreTest, Count) {
    // assume that ds.set() and add() work fine
    EXPECT_EQ(8, ds0_->count());

    // If the DataStore is not inisialized yet, it returns 0.
    kmrnext::DataStore *ds1 = gNext->create_ds(ds_size_);
    EXPECT_EQ(0, ds1->count());
    ds1->set(array_ds0_);
    EXPECT_EQ(0, ds1->count());

    // If one Data is added, it returns 1.
    ds1->add(*key0_, *d0_);
    EXPECT_EQ(1, ds1->count());
    delete ds1;
  }

  TEST_F(DataStoreTest, Zeroize) {
    kmrnext::DataStore *ds = gNext->create_ds(ds_size_);
    ds->set(array_ds0_);
    ds->zeroize();
    kmrnext::Key key(ds_size_);
    for (size_t i = 0; i < array_ds0_[0]; i++) {
      for (size_t j = 0; j < array_ds0_[1]; j++) {
	for (size_t k = 0; k < array_ds0_[2]; k++) {
	  key.set_dim(0, i);
	  key.set_dim(1, j);
	  key.set_dim(2, k);
	  kmrnext::DataPack dp = ds->get(key);
	  EXPECT_EQ(0, *static_cast<long*>(dp.data().value()));
	  EXPECT_EQ(sizeof(long), dp.data().size());
	}
      }
    }
    delete ds;
  }

  TEST_F(DataStoreTest, Duplicate) {
    kmrnext::DataStore* ds0 = ds0_->duplicate();
    // Pointer to the duplicated DataStore is different
    EXPECT_NE(ds0_, ds0);
    // Pointers to data elements having the same coordinate are not same
    kmrnext::DataPack dp0_ = ds0_->get(*key0_);
    kmrnext::DataPack dp0  = ds0 ->get(*key0_);
    void *v_ds0_ = dp0_.data().value();
    void *v_ds0  = dp0 .data().value();
    // The following codes do not work as I thought, why?
    // Maybe, memory of v_ds0_ is freed when the second get() is called.
    //   void *v_ds0_ = ds0_->get(*key0_).data()->value();
    //   void *v_ds0  = ds0 ->get(*key0_).data()->value();
    EXPECT_NE(v_ds0_, v_ds0);
    // Values of data elements having the same coordinate are same
    EXPECT_EQ(*static_cast<long*>(v_ds0_), *static_cast<long*>(v_ds0));
    delete ds0;
  }

}
