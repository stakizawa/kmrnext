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
	ds2_vec.push_back((int)i + 1);
      }
      ds2_owners_ = new int[4];
      init_owners(ds2_owners_, 4);
      DataLoader1D ds2_loader1d(ds2_array_[1]);
      ds2_ = new kmrnext::DataStore(2, gNext);
      ds2_->set(ds2_array_);
      ds2_->load_array(ds2_vec, ds2_loader1d);

      ds3_array_ = new size_t[3];
      ds3_array_[0] = 4;
      ds3_array_[1] = 4;
      ds3_array_[2] = 4;
      std::vector<int> ds3_vec;
      for (size_t i = 0; i < ds3_array_[0] * ds3_array_[1]; i++) {
	ds3_vec.push_back(1);
      }
      DataLoader1D ds3_loader1d(ds3_array_[2]);
      ds3_ = new kmrnext::DataStore(3, gNext);
      ds3_->set(ds3_array_);
      ds3_->load_array(ds3_vec, ds3_loader1d);


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
    // If the coordinates whose View is Ture are same,
    // the owner process is same.
    kmrnext::DataStore ods0(3, gNext);
    ods0.set(ds3_array_);
    SummarizerSingle0 mapper0(rank);
    kmrnext::View v0(3);
    bool flags0[3] = {true, false, true};
    v0.set(flags0);
    ds3_->map_single(&ods0, mapper0, v0);
    EXPECT_EQ(ods0.get(*k3_000_).data()->owner(),
	      ods0.get(*k3_030_).data()->owner());

    kmrnext::DataStore ods1(3, gNext);
    ods1.set(ds3_array_);
    kmrnext::View v1(3);
    bool flags1[3] = {false, true, false};
    v1.set(flags1);
    ds3_->map_single(&ods1, mapper0, v1);
    EXPECT_EQ(ods1.get(*k3_010_).data()->owner(),
	      ods1.get(*k3_113_).data()->owner());

    // If all fields of a View is false, all data is gathered to a
    // specific rank
    kmrnext::DataStore ods2(3, gNext);
    ods2.set(ds3_array_);
    kmrnext::View v2(3);
    bool flags2[3] = {false, false, false};
    v2.set(flags2);
    ds3_->map_single(&ods2, mapper0, v2);
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
    ds3_->map_single(&ods3, mapper0, v3);
  }

}
