// This test tests DataStore class for KMR backend.
#include <gtest/gtest.h>
#include <sstream>
#include <cassert>
#include <mpi.h>
#include "kmrnext.hpp"

extern kmrnext::KMRNext *gNext;

namespace {

  class DataLoader2D : public kmrnext::DataStore::Loader<long> {
    size_t size_;
  public:
    DataLoader2D(size_t siz) : size_(siz) {}

    int operator()(kmrnext::DataStore *ds, const long& num)
    {
      kmrnext::Key key(2);
      key.set_dim(0, num);
      int val = (int)num + 1;
      kmrnext::Data data((void*)&val, sizeof(int));
      for (size_t i = 0; i < size_; i++) {
	key.set_dim(1, i);
	ds->add(key, data);
      }
      return 0;
    }
  };

  class DataLoader3D : public kmrnext::DataStore::Loader<long> {
    size_t size_y_;
    size_t size_z_;
  public:
    DataLoader3D(size_t siz_y, size_t siz_z)
      : size_y_(siz_y), size_z_(siz_z) {}

    int operator()(kmrnext::DataStore *ds, const long& num)
    {
      size_t x = num / size_y_;
      size_t y = num % size_y_;
      kmrnext::Key key(3);
      key.set_dim(0, x);
      key.set_dim(1, y);
      int val = 1;
      kmrnext::Data data((void*)&val, sizeof(int));
      for (size_t i = 0; i < size_z_; i++) {
	key.set_dim(2, i);
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
      std::vector<long> ds2_vec;
      for (size_t i = 0; i < ds2_array_[0]; i++) {
	ds2_vec.push_back(i);
      }
      ds2_owners_ = new int[4];
      init_owners(ds2_owners_, 4);
      DataLoader2D ds2_loader(ds2_array_[1]);
      ds2_ = new kmrnext::DataStore(2, gNext);
      ds2_->set(ds2_array_);
      ds2_->load_integers(ds2_vec, ds2_loader);

      ds3_array_ = new size_t[3];
      ds3_array_[0] = 4;
      ds3_array_[1] = 4;
      ds3_array_[2] = 4;
      std::vector<long> ds3_vec;
      for (size_t i = 0; i < ds3_array_[0] * ds3_array_[1]; i++) {
	ds3_vec.push_back(i);
      }
      DataLoader3D ds3_loader(ds3_array_[1], ds3_array_[2]);
      ds3_ = new kmrnext::DataStore(3, gNext);
      ds3_->set(ds3_array_);
      ds3_->load_integers(ds3_vec, ds3_loader);


      size_t ary_k2_00[2] = {0,0};
      init_key(&k2_00_, 2, ary_k2_00);
      size_t ary_k2_11[2] = {1,1};
      init_key(&k2_11_, 2, ary_k2_11);
      size_t ary_k2_22[2] = {2,2};
      init_key(&k2_22_, 2, ary_k2_22);
      size_t ary_k2_33[2] = {3,3};
      init_key(&k2_33_, 2, ary_k2_33);

      size_t ary_k3_000[3] = {0,0,0};
      init_key(&k3_000_, 3, ary_k3_000);
      size_t ary_k3_030[3] = {0,3,0};
      init_key(&k3_030_, 3, ary_k3_030);
      size_t ary_k3_010[3] = {0,1,0};
      init_key(&k3_010_, 3, ary_k3_010);
      size_t ary_k3_113[3] = {1,1,3};
      init_key(&k3_113_, 3, ary_k3_113);
    }

    void init_key(kmrnext::Key **kp, size_t size, size_t *array) {
      *kp = new kmrnext::Key(size);
      (*kp)->set(array);
    }

    void init_owners(int *owners, size_t size) {
      int quotient = (int)(size / (size_t)nprocs);
      int remain   = (int)(size % (size_t)nprocs);
      for (int i = 0; i < nprocs; i++) {
	int start = i * quotient + ((i < remain)? i : remain);
	int end = start + quotient + ((i < remain)? 1 : 0);
	if ((size_t)start >= size) { break; }
	for (int j = start; j < end; j++) {
	  owners[j] = i;
	}
      }
    }

    // Create a DataStore <3,3,3> having '1' for all data element
    kmrnext::DataStore *create_ds3d() {
      kmrnext::DataStore *ds3 = new kmrnext::DataStore(3, gNext);
      ds3->set_dim(0, 3);
      ds3->set_dim(1, 3);
      ds3->set_dim(2, 3);
      kmrnext::Key ds3_key(3);
      int ds3_val = 1;
      kmrnext::Data ds3_dat(&ds3_val, sizeof(int));
      for (size_t i = 0; i < 3; i++) {
	ds3_key.set_dim(0, i);
	for (size_t j = 0; j < 3; j++) {
	  ds3_key.set_dim(1, j);
	  for (size_t k = 0; k < 3; k++) {
	    ds3_key.set_dim(2, k);
	    ds3->add(ds3_key, ds3_dat);
	  }
	}
      }
      return ds3;
    }

    virtual ~KMRDataStoreTest() {
      delete[] ds2_array_;
      delete   ds2_;
      delete[] ds2_owners_;
      delete[] ds3_array_;
      delete   ds3_;
      delete   k2_00_;
      delete   k2_11_;
      delete   k2_22_;
      delete   k2_33_;
      delete   k3_000_;
      delete   k3_030_;
      delete   k3_010_;
      delete   k3_113_;
    }

    int rank;
    int nprocs;

    size_t *ds2_array_;       // {4,4}
    kmrnext::DataStore *ds2_; // Dimension: {4,4}
                              // Data: <0,x> => 1, <1,x> => 2, <2,x> => 3
                              //       <3,x> => 4
    int *ds2_owners_;         // size is 4, calculated at runtime

    size_t *ds3_array_;       // {4,4,4}
    kmrnext::DataStore *ds3_; // Dimension: {4,4,4}
                              // Data: all 1

    kmrnext::Key *k2_00_;     // <0,0>
    kmrnext::Key *k2_11_;     // <0,0>
    kmrnext::Key *k2_22_;     // <2,2>
    kmrnext::Key *k2_33_;     // <3,3>

    kmrnext::Key *k3_000_;    // <0,0,0>
    kmrnext::Key *k3_030_;    // <0,3,0>
    kmrnext::Key *k3_010_;    // <0,1,0>
    kmrnext::Key *k3_113_;    // <1,1,3>
  };

  TEST_F(KMRDataStoreTest, Load_integers) {
    size_t ds_size[2] = {4,4};
    int owners[4];
    DataLoader2D loader(ds_size[1]);

    std::vector<long> vec0;
    for (size_t i = 0; i < ds_size[0]; i++) {
      vec0.push_back(i);
    }
    init_owners(owners, 4);
    kmrnext::DataStore ds0(2, gNext);
    ds0.set(ds_size);
    ds0.load_integers(vec0, loader);
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
    ds2_->map(mapper0, v0, &ods0);
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
    ds2_->map(mapper1, v1, &ods1);
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

  TEST_F(KMRDataStoreTest, Set_allocation_view) {
    kmrnext::View ds3_0_dav(3);
    bool flags3_0_dav[3] = {true, false, false};
    ds3_0_dav.set(flags3_0_dav);

    // If a DataStore is initialized without calling load_xxx(),
    // the Allocation View is set <T, F, ...>, by default.
    kmrnext::DataStore ds3_0(3, gNext);
    {
      ds3_0.set_dim(0, 2);
      ds3_0.set_dim(1, 2);
      ds3_0.set_dim(2, 2);
      kmrnext::Key ds3_0_key(3);
      int ds3_0_val = 1;
      kmrnext::Data ds3_0_dat(&ds3_0_val, sizeof(int));
      for (size_t i = 0; i < 2; i++) {
	ds3_0_key.set_dim(0, i);
	for (size_t j = 0; j < 2; j++) {
	  ds3_0_key.set_dim(1, j);
	  for (size_t k = 0; k < 2; k++) {
	    ds3_0_key.set_dim(2, k);
	    ds3_0.add(ds3_0_key, ds3_0_dat);
	  }
	}
      }
    }
    EXPECT_EQ(ds3_0_dav, ds3_0.get_allocation_view());

    kmrnext::View v3_0(3);
    bool flags3_0[3] = {true, true, false};
    v3_0.set(flags3_0);
    kmrnext::View v3_1(3);
    bool flags3_1[3] = {false, true, true};
    v3_1.set(flags3_1);
    kmrnext::View v3_2(3);
    bool flags3_2[3] = {false, false, true};
    v3_2.set(flags3_2);
    kmrnext::View v2(2);
    bool flags2[2] = {true, true};
    v2.set(flags2);

    // If a DataStore is initialized by calling load_xxx(),
    // the Allocation View is set according to the number of array items.
    kmrnext::View ds3_dav = ds3_->get_allocation_view();
    EXPECT_EQ(v3_0, ds3_dav);

    // If a view is set to the DataStore, the gotten View should be same.
    ds3_->set_allocation_view(v3_1);
    kmrnext::View v0 = ds3_->get_allocation_view();
    EXPECT_EQ(v3_1, v0);

    // If another view is set to the DataStore, the View should be replaced.
    ds3_->set_allocation_view(v3_2);
    v0 = ds3_->get_allocation_view();
    EXPECT_EQ(v3_2, v0);

    // If the size of DataStore and View is not same, it throws runtime_error.
    EXPECT_THROW({ds3_->set_allocation_view(v2);}, std::runtime_error);
  }

  TEST_F(KMRDataStoreTest, Collate) {
    kmrnext::Key key000(3);
    size_t _k000[3] = {0, 0, 0};
    key000.set(_k000);
    kmrnext::Key key111(3);
    size_t _k111[3] = {1, 1, 1};
    key111.set(_k111);
    kmrnext::Key key222(3);
    size_t _k222[3] = {2, 2, 2};
    key222.set(_k222);
    kmrnext::Key key012(3);
    size_t _k012[3] = {0, 1, 2};
    key012.set(_k012);
    kmrnext::Key key210(3);
    size_t _k210[3] = {2, 1, 0};
    key210.set(_k210);

    kmrnext::DataStore *ds0 = create_ds3d();
    // As the added data without calling map() reside on rank 0 only, and
    // the default Allocation View is <T, F, F>, the first calling collate()
    // makes change of data allocation.
    // Check the allocation before collate()
    EXPECT_EQ(0, ds0->get(key000).data()->owner());
    EXPECT_EQ(0, ds0->get(key111).data()->owner());
    EXPECT_EQ(0, ds0->get(key222).data()->owner());
    ds0->collate();
    // Check the allocation after collate()
    if (nprocs >= 3) {
      EXPECT_EQ(0, ds0->get(key000).data()->owner());
      EXPECT_EQ(1, ds0->get(key111).data()->owner());
      EXPECT_EQ(2, ds0->get(key222).data()->owner());
    } else if (nprocs == 2) {
      EXPECT_EQ(0, ds0->get(key000).data()->owner());
      EXPECT_EQ(0, ds0->get(key111).data()->owner());
      EXPECT_EQ(1, ds0->get(key222).data()->owner());
    } else {  // nprocs == 1
      EXPECT_EQ(0, ds0->get(key000).data()->owner());
      EXPECT_EQ(0, ds0->get(key111).data()->owner());
      EXPECT_EQ(0, ds0->get(key222).data()->owner());
    }
    EXPECT_EQ(1, *(int*)ds0->get(key000).data()->value());
    EXPECT_EQ(1, *(int*)ds0->get(key111).data()->value());
    EXPECT_EQ(1, *(int*)ds0->get(key222).data()->value());

    // Without resetting the Allocation View, calling collate() again does
    // not take any effect.
    ds0->collate();
    if (nprocs >= 3) {
      EXPECT_EQ(0, ds0->get(key000).data()->owner());
      EXPECT_EQ(1, ds0->get(key111).data()->owner());
      EXPECT_EQ(2, ds0->get(key222).data()->owner());
    } else if (nprocs == 2) {
      EXPECT_EQ(0, ds0->get(key000).data()->owner());
      EXPECT_EQ(0, ds0->get(key111).data()->owner());
      EXPECT_EQ(1, ds0->get(key222).data()->owner());
    } else {  // nprocs == 1
      EXPECT_EQ(0, ds0->get(key000).data()->owner());
      EXPECT_EQ(0, ds0->get(key111).data()->owner());
      EXPECT_EQ(0, ds0->get(key222).data()->owner());
    }

    kmrnext::View aviewFTT(3);
    bool _avFTT[3] = {false, true, true};
    aviewFTT.set(_avFTT);

    // By changing the Allocation View, calling collate() makes change of
    // data allocation.
    // The resultant viewed keys by applying this Allocation View is listed
    // below and each of them has 3 data.
    //
    //        Keys  nprocs 1  2  3  4  5  6  7  8  9
    //   <*, 0, 0>         0  0  0  0  0  0  0  0  0
    //   <*, 0, 1>         0  0  0  0  0  0  0  0  1
    //   <*, 0, 2>         0  0  0  0  1  1  1  1  2
    //   <*, 1, 0>         0  0  1  1  1  1  1  2  3
    //   <*, 1, 1>         0  0  1  1  2  2  2  3  4
    //   <*, 1, 2>         0  1  1  2  2  2  3  4  5
    //   <*, 2, 0>         0  1  2  2  3  3  4  5  6
    //   <*, 2, 1>         0  1  2  3  3  4  5  6  7
    //   <*, 2, 2>         0  1  2  3  4  5  6  7  8
    ds0->set_allocation_view(aviewFTT);
    ds0->collate();
    if (nprocs >= 9) {
      EXPECT_EQ(0, ds0->get(key000).data()->owner());
      EXPECT_EQ(4, ds0->get(key111).data()->owner());
      EXPECT_EQ(8, ds0->get(key222).data()->owner());
      EXPECT_EQ(5, ds0->get(key012).data()->owner());
      EXPECT_EQ(3, ds0->get(key210).data()->owner());
    } else if (nprocs == 8) {
      EXPECT_EQ(0, ds0->get(key000).data()->owner());
      EXPECT_EQ(3, ds0->get(key111).data()->owner());
      EXPECT_EQ(7, ds0->get(key222).data()->owner());
      EXPECT_EQ(4, ds0->get(key012).data()->owner());
      EXPECT_EQ(2, ds0->get(key210).data()->owner());
    } else if (nprocs == 7) {
      EXPECT_EQ(0, ds0->get(key000).data()->owner());
      EXPECT_EQ(2, ds0->get(key111).data()->owner());
      EXPECT_EQ(6, ds0->get(key222).data()->owner());
      EXPECT_EQ(3, ds0->get(key012).data()->owner());
      EXPECT_EQ(1, ds0->get(key210).data()->owner());
    } else if (nprocs == 6) {
      EXPECT_EQ(0, ds0->get(key000).data()->owner());
      EXPECT_EQ(2, ds0->get(key111).data()->owner());
      EXPECT_EQ(5, ds0->get(key222).data()->owner());
      EXPECT_EQ(2, ds0->get(key012).data()->owner());
      EXPECT_EQ(1, ds0->get(key210).data()->owner());
    } else if (nprocs == 5) {
      EXPECT_EQ(0, ds0->get(key000).data()->owner());
      EXPECT_EQ(2, ds0->get(key111).data()->owner());
      EXPECT_EQ(4, ds0->get(key222).data()->owner());
      EXPECT_EQ(2, ds0->get(key012).data()->owner());
      EXPECT_EQ(1, ds0->get(key210).data()->owner());
    } else if (nprocs == 4) {
      EXPECT_EQ(0, ds0->get(key000).data()->owner());
      EXPECT_EQ(1, ds0->get(key111).data()->owner());
      EXPECT_EQ(3, ds0->get(key222).data()->owner());
      EXPECT_EQ(2, ds0->get(key012).data()->owner());
      EXPECT_EQ(1, ds0->get(key210).data()->owner());
    } else if (nprocs == 3) {
      EXPECT_EQ(0, ds0->get(key000).data()->owner());
      EXPECT_EQ(1, ds0->get(key111).data()->owner());
      EXPECT_EQ(2, ds0->get(key222).data()->owner());
      EXPECT_EQ(1, ds0->get(key012).data()->owner());
      EXPECT_EQ(1, ds0->get(key210).data()->owner());
    } else if (nprocs == 2) {
      EXPECT_EQ(0, ds0->get(key000).data()->owner());
      EXPECT_EQ(0, ds0->get(key111).data()->owner());
      EXPECT_EQ(1, ds0->get(key222).data()->owner());
      EXPECT_EQ(1, ds0->get(key012).data()->owner());
      EXPECT_EQ(0, ds0->get(key210).data()->owner());
    } else {  // nprocs == 1
      EXPECT_EQ(0, ds0->get(key000).data()->owner());
      EXPECT_EQ(0, ds0->get(key111).data()->owner());
      EXPECT_EQ(0, ds0->get(key222).data()->owner());
      EXPECT_EQ(0, ds0->get(key012).data()->owner());
      EXPECT_EQ(0, ds0->get(key210).data()->owner());
    }
  }

}
