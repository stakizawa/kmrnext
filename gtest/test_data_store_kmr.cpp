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

#if 0
  TEST_F(KMRDataStoreTest, Load_array) {
    // also tests ds.get()
    size_t ds_size[2] = {4,4};
    int owners[4];
    DataLoader1D loader1d(ds_size[1]);

    // Correctly load data by each process
    std::vector<int> vec0;
    for (size_t i = 0; i < ds_size[0]; i++) {
      vec0.push_back((int)i + 1);
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

  // A mapper class that calculates sum of all data.  The resultant value
  // is written to the coordinates where summed values are located.
  class SummarizerSingle0 : public kmrnext::DataStore::Mapper {
    int base_;
  public:
    SummarizerSingle0(int base) : base_((base + 1) * 1000) {}
    int operator()(kmrnext::DataStore *inds, kmrnext::DataStore *outds,
		   kmrnext::Key& key, std::vector<kmrnext::DataPack>& dps,
		   kmrnext::DataStore::MapEnvironment& env)
    {
      {
	int nprocs;
	MPI_Comm_size(env.mpi_comm, &nprocs);
	assert(nprocs == 1);
      }
      int val = base_;
      for (std::vector<kmrnext::DataPack>:: iterator itr = dps.begin();
	   itr != dps.end(); itr++) {
	val += *(int*)(*itr).data()->value();
      }
      kmrnext::Data d(&val, sizeof(int));

      for (std::vector<kmrnext::DataPack>:: iterator itr = dps.begin();
	   itr != dps.end(); itr++) {
	kmrnext::Key k = (*itr).key();
	outds->add(k, d);
      }
      return 0;
    }
  };

  TEST_F(KMRDataStoreTest, Map_single) {
    // The owners of data in the output DS whose coordinates are
    // same by applying the View are same.
    kmrnext::DataStore ods0(3, gNext);
    ods0.set(ds3_array_);
    SummarizerSingle0 mapper0(rank);
    kmrnext::View v0(3);
    bool flags0[3] = {true, false, true};
    v0.set(flags0);
    ds3_->map_single(mapper0, v0, &ods0);
    EXPECT_EQ(ods0.get(*k3_000_).data()->owner(),
	      ods0.get(*k3_030_).data()->owner());

    kmrnext::DataStore ods1(3, gNext);
    ods1.set(ds3_array_);
    kmrnext::View v1(3);
    bool flags1[3] = {false, true, false};
    v1.set(flags1);
    ds3_->map_single(mapper0, v1, &ods1);
    EXPECT_EQ(ods1.get(*k3_010_).data()->owner(),
	      ods1.get(*k3_113_).data()->owner());

    // If all fields of a View is false, all data is gathered to a
    // specific rank
    kmrnext::DataStore ods2(3, gNext);
    ods2.set(ds3_array_);
    kmrnext::View v2(3);
    bool flags2[3] = {false, false, false};
    v2.set(flags2);
    ds3_->map_single(mapper0, v2, &ods2);
    EXPECT_EQ(ods2.get(*k3_000_).data()->owner(),
	      ods2.get(*k3_030_).data()->owner());
    EXPECT_EQ(ods2.get(*k3_010_).data()->owner(),
	      ods2.get(*k3_113_).data()->owner());
    EXPECT_EQ(ods2.get(*k3_000_).data()->owner(),
	      ods2.get(*k3_010_).data()->owner());

    // If all field of a View is true, DataStore::map() is called.
    kmrnext::DataStore ods3(3, gNext);
    ods3.set(ds3_array_);
    kmrnext::View v3(3);
    bool flags3[3] = {true, true, true};
    v3.set(flags3);
    ds3_->map_single(mapper0, v3, &ods3);
    EXPECT_EQ(ds3_->get(*k3_000_).data()->owner(),
	      ods3.get(*k3_000_).data()->owner());
    EXPECT_EQ(ds3_->get(*k3_030_).data()->owner(),
	      ods3.get(*k3_030_).data()->owner());
    EXPECT_EQ(ds3_->get(*k3_113_).data()->owner(),
	      ods3.get(*k3_113_).data()->owner());

    // If the output DS is omitted, the results are writted to the input DS
    kmrnext::DataStore *ds3 = ds3_->duplicate();
    ds3->map_single(mapper0, v0);
    EXPECT_EQ(*(int*)ods0.get(*k3_000_).data()->value(),
	      *(int*)ds3->get(*k3_000_).data()->value());
    EXPECT_EQ(*(int*)ods0.get(*k3_030_).data()->value(),
	      *(int*)ds3->get(*k3_030_).data()->value());
    EXPECT_EQ(ds3->get(*k3_000_).data()->owner(),
	      ds3->get(*k3_030_).data()->owner());
    delete ds3;

    ds3 = ds3_->duplicate();
    ds3->map_single(mapper0, v1);
    EXPECT_EQ(*(int*)ods1.get(*k3_010_).data()->value(),
    	      *(int*)ds3->get(*k3_010_).data()->value());
    EXPECT_EQ(ds3->get(*k3_010_).data()->owner(),
     	      ds3->get(*k3_113_).data()->owner());
    delete ds3;

    ds3 = ds3_->duplicate();
    ds3->map_single(mapper0, v2);
    EXPECT_EQ(*(int*)ods2.get(*k3_000_).data()->value(),
    	      *(int*)ds3->get(*k3_000_).data()->value());
    EXPECT_EQ(*(int*)ods2.get(*k3_113_).data()->value(),
    	      *(int*)ds3->get(*k3_113_).data()->value());
    EXPECT_EQ(ds3->get(*k3_000_).data()->owner(),
    	      ds3->get(*k3_030_).data()->owner());
    EXPECT_EQ(ds3->get(*k3_010_).data()->owner(),
    	      ds3->get(*k3_113_).data()->owner());
    delete ds3;
  }

  TEST_F(KMRDataStoreTest, Collate) {
    if (nprocs < 4) {
      EXPECT_TRUE(false) << "Test for DataStore.collate() is skipped "
			 << "as there are not enough number of processes.  "
			 << "Specify more than 4 processes.";
      return;
    }

    // The same operation of the first test in Map_single.
    kmrnext::DataStore ods0(3, gNext);
    ods0.set(ds3_array_);
    SummarizerSingle0 mapper0(rank);
    kmrnext::View v0(3);
    bool flags0[3] = {true, false, true};
    v0.set(flags0);
    ds3_->map_single(mapper0, v0, &ods0);
    EXPECT_EQ(ods0.get(*k3_000_).data()->owner(),
	      ods0.get(*k3_030_).data()->owner());

    // 1. Data whose coorinate of the second dimension are same are
    //    gathered to the same node.
    kmrnext::View cv0(3);
    bool cflags0[3] = {false, true, false};
    cv0.set(cflags0);
    ods0.collate(cv0);
    EXPECT_EQ(ods0.get(*k3_010_).data()->owner(),
	      ods0.get(*k3_113_).data()->owner());
    EXPECT_NE(ods0.get(*k3_000_).data()->owner(),
	      ods0.get(*k3_030_).data()->owner());

    // 2. Data whose coorinate of the first dimension are same aregathered
    //    gathered to the same node.
    kmrnext::View cv1(3);
    bool cflags1[3] = {true, false, false};
    cv1.set(cflags1);
    ods0.collate(cv1);
    EXPECT_EQ(ods0.get(*k3_010_).data()->owner(),
	      ods0.get(*k3_000_).data()->owner());

    // The ame operation of the first test in Map_single, except that
    // it performs in-place mapping.
    kmrnext::DataStore *ds3 = ds3_->duplicate();
    ds3->map_single(mapper0, v0);
    EXPECT_EQ(*(int*)ods0.get(*k3_000_).data()->value(),
	      *(int*)ds3->get(*k3_000_).data()->value());
    EXPECT_EQ(*(int*)ods0.get(*k3_030_).data()->value(),
	      *(int*)ds3->get(*k3_030_).data()->value());
    EXPECT_EQ(ds3->get(*k3_000_).data()->owner(),
	      ds3->get(*k3_030_).data()->owner());

    // 3. Data whose coorinate of the second dimension are same are
    //    gathered to the same node.
    ds3->collate(cv0);
    EXPECT_EQ(ds3->get(*k3_010_).data()->owner(),
	      ds3->get(*k3_113_).data()->owner());
    EXPECT_NE(ds3->get(*k3_000_).data()->owner(),
	      ds3->get(*k3_030_).data()->owner());

    // 4. Data whose coorinate of the first dimension are same aregathered
    //    gathered to the same node.
    ds3->collate(cv1);
    EXPECT_EQ(ds3->get(*k3_010_).data()->owner(),
	      ds3->get(*k3_000_).data()->owner());
    delete ds3;
  }
#endif

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

}
