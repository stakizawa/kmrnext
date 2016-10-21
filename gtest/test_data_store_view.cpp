// This test tests DataStore functions using detailed views.
// The Detailed view is a view that has numbers in it to specify the sizes of
// blocks that split each dimension.
// Some tests in this test is only available in KMR backend.
#include <gtest/gtest.h>
#ifdef BACKEND_KMR
#include <mpi.h>
#endif
#include "kmrnext.hpp"

extern kmrnext::KMRNext *gNext;

namespace {

  class DataStoreViewTest : public ::testing::Test {
  protected:
    DataStoreViewTest() {
#ifdef BACKEND_KMR
      MPI_Comm_rank(gNext->kmr()->comm, &rank);
      MPI_Comm_size(gNext->kmr()->comm, &nprocs);
#else
      rank   = 0;
      nprocs = 1;
#endif

      // Create a 2D DataStore.
      // In case of KMR backend, the following procedure create an
      // instance of 2D DataStore on each process, but the value(1)
      // are put only on rank 0 process.
      ds_ = gNext->create_ds(2);
      size_t dim_siz = nprocs * 2;
      ds_->set_dim(0, dim_siz);
      ds_->set_dim(1, dim_siz);
      kmrnext::Key key(2);
      int value = 1;
      kmrnext::Data dat(&value, sizeof(int));
      for (size_t i = 0; i < dim_siz; i++) {
	key.set_dim(0, i);
	for (size_t j = 0; j < dim_siz; j++) {
	  key.set_dim(1, j);
	  ds_->add(key, dat);
	}
      }

      view_ = create_view();
      split_ = create_view();
    }

    virtual ~DataStoreViewTest() {
      delete ds_;
      delete view_;
      delete split_;
    }

    // Create a view/split
    kmrnext::View* create_view() {
      long ary[2];
      if (nprocs > 2 && nprocs % 2 == 0) {
	// Divide the second dimension by 2 and the first dimension
	// by using remaining processes
	ary[0] = nprocs / 2;
	ary[1] = 2;
      } else {
	// Do not divie the first dimension, but divide the second
	// dimension by all processes
	ary[0] = kmrnext::View::SplitNone;
	ary[1] = nprocs;
      }
      kmrnext::View* v = new kmrnext::View(2);
      v->set(ary);
      return v;
    }

    // Returns count of data elements in a group.
    // It counts NULL data elements.
    size_t count_in_block(kmrnext::View* v) {
      size_t count;
      size_t dim_siz = 2 * nprocs;
      if (nprocs > 2 && nprocs % 2 == 0) {
	count = (dim_siz / v->dim(0)) * (dim_siz / v->dim(1));
      } else {
	count = dim_siz * (dim_siz / v->dim(1));
      }
      return count;
    }

    // Create an onwer list that has owners of each data elements defined
    // by a Split created by create_view()
    int* create_owners() {
      size_t dim_siz = nprocs * 2;
      int* owners = new int[dim_siz * dim_siz];

      // array is same as the above
      long blk[2];
      if (nprocs > 2 && nprocs % 2 == 0) {
	blk[0] = 4;       // 2 * nprocs / ary[0]
	blk[1] = nprocs;  // 2 * nprocs / 2
#if 1
	// The following is the case of column-ordered index
	size_t block_i = dim_siz / blk[0];
	for (size_t i = 0; i < dim_siz; i++) {
	  size_t chnk_i = i / blk[0];
	  for (size_t j = 0; j < dim_siz; j++) {
	    size_t chnk_j = j / blk[1];
	    int owner = static_cast<int>(chnk_j * block_i + chnk_i);
	    owners[j * dim_siz + i] = owner;
	  }
	}
#else
	// The following is the case of row-ordered index
	size_t block_j = dim_siz / blk[1];
	for (size_t i = 0; i < dim_siz; i++) {
	  size_t chnk_i = i / blk[0];
	  for (size_t j = 0; j < dim_siz; j++) {
	    size_t chnk_j = j / blk[1];
	    int owner = static_cast<int>(chnk_i * block_j + chnk_j);
	    owners[j * dim_siz + i] = owner;
	  }
	}
#endif
      } else {
	blk[0] = 2 * nprocs;  // as do not split
	blk[1] = 2;           // n * nprocs / ary[1]
	for (size_t i = 0; i < dim_siz; i++) {
	  for (size_t j = 0; j < dim_siz; j++) {
	    int owner = static_cast<int>(j / blk[1]);
	    owners[j * dim_siz + i] = owner;
	  }
	}
      }

      return owners;
    }

    void delete_owners(int *owners) {
      delete[] owners;
    }

    int rank;
    int nprocs;

    kmrnext::DataStore* ds_;  // Dimension:  {nprocs*2, nprocs*2}
    kmrnext::View* view_;     // If number of nprocs is even, <nprocs/2, 2>
                              // If number of nprocs is odd, <NONE, nprocs>
    kmrnext::View* split_;    // same as view_
  };

  TEST_F(DataStoreViewTest, Basic) {
    // A runtime error is thrown when the split count of at least one
    // dimension in a View/Split is too big and can not divide the
    // size of dimensions in a DataStore

    kmrnext::View view0(2);
    long ary_too_big[2] = {nprocs * 100, nprocs * 100};
    view0.set(ary_too_big);
    kmrnext::Key key0(2);
    size_t ary_key[2] = {0, 0};
    key0.set(ary_key);
    EXPECT_THROW({ds_->get(view0, key0);}, std::runtime_error);

#ifdef BACKEND_KMR
    kmrnext::View split0(2);
    split0.set(ary_too_big);
    EXPECT_THROW({ds_->set_split(split0);}, std::runtime_error);
#endif
  }

  TEST_F(DataStoreViewTest, GetView) {
    // length of a block
    size_t len0, len1;
    {
      if (nprocs > 2 && nprocs % 2 == 0) {
	len0 = ds_->dim(0) / view_->dim(0);
      } else {
	len0 = ds_->dim(0);
      }
      len1 = ds_->dim(1) / view_->dim(1);
    }

    // Get data elements with key<0, 0> and the view
    kmrnext::Key k00(2);
    k00.set_dim(0, 0);
    k00.set_dim(1, 0);
    size_t k00_0_head = k00.dim(0);
    size_t k00_0_tail = k00.dim(0) + len0;
    size_t k00_1_head = k00.dim(1);
    size_t k00_1_tail = k00.dim(1) + len1;
    std::vector<kmrnext::DataPack> *vec0 = ds_->get(*view_, k00);
    // Check if the gotten data count is correct
    EXPECT_EQ(count_in_block(view_), vec0->size());
    // Check if keys of gotten data is in the correct region
    for (std::vector<kmrnext::DataPack>::iterator itr = vec0->begin();
	 itr != vec0->end(); itr++) {
      kmrnext::Key k = (*itr).key();
      EXPECT_TRUE(k00_0_head <= k.dim(0) && k.dim(0) < k00_0_tail);
      EXPECT_TRUE(k00_1_head <= k.dim(1) && k.dim(1) < k00_1_tail);
    }
    delete vec0;

    // Get data elements with key<LAST, LAST> and the view
    kmrnext::Key kll(2);
    kll.set_dim(0, ds_->dim(0) / len0 - 1);
    kll.set_dim(1, ds_->dim(1) / len1 - 1);
    size_t kll_0_head = kll.dim(0) * len0;
    size_t kll_0_tail = kll_0_head + len0;
    size_t kll_1_head = kll.dim(1) * len1;
    size_t kll_1_tail = kll_1_head + len1;
    std::vector<kmrnext::DataPack> *vec1 = ds_->get(*view_, kll);
    // Check if the gotten data count is correct
    EXPECT_EQ(count_in_block(view_), vec1->size());
    // Check if keys of gotten data is in the correct region
    for (std::vector<kmrnext::DataPack>::iterator itr = vec1->begin();
	 itr != vec1->end(); itr++) {
      kmrnext::Key k = (*itr).key();
      EXPECT_TRUE(kll_0_head <= k.dim(0) && k.dim(0) < kll_0_tail);
      EXPECT_TRUE(kll_1_head <= k.dim(1) && k.dim(1) < kll_1_tail);
    }
    delete vec1;
  }

  class MapTestMapper0 : public kmrnext::DataStore::Mapper {
    int nprocs_;
  public:
    MapTestMapper0(int nprocs) : nprocs_(nprocs) {}

    int operator()(kmrnext::DataStore *inds, kmrnext::DataStore *outds,
		   kmrnext::Key& key, std::vector<kmrnext::DataPack>& dps,
		   kmrnext::DataStore::MapEnvironment& env)
    {
      size_t blk_len0, blk_len1;
      {
	if (nprocs_ > 2 && nprocs_ % 2 == 0) {
	  blk_len0 = inds->dim(0) / env.view.dim(0);
	} else {
	  blk_len0 = inds->dim(0);
	}
	blk_len1 = inds->dim(1) / env.view.dim(1);
      }

#ifdef BACKEND_SERIAL
      EXPECT_EQ(blk_len0 * blk_len1, dps.size());
#endif
#ifdef BACKEND_KMR
      int local_count = static_cast<int>(dps.size());
      int total_count;
      MPI_Allreduce(&local_count, &total_count, 1, MPI_INT, MPI_SUM,
		    env.mpi_comm);
      EXPECT_EQ(blk_len0 * blk_len1, total_count);
#endif

      size_t k_0_head, k_0_tail, k_1_head, k_1_tail;
      {
	kmrnext::DataPack dp0 = dps.at(0);
	kmrnext::Key k = dp0.key();
	k_0_head = (k.dim(0) / blk_len0) * blk_len0;
	k_0_tail = k_0_head + blk_len0;
	k_1_head = (k.dim(1) / blk_len1) * blk_len1;
	k_1_tail = k_1_head + blk_len1;
      }
      for (std::vector<kmrnext::DataPack>::iterator itr = dps.begin();
	   itr != dps.end(); itr++) {
	kmrnext::Key k = (*itr).key();
	EXPECT_TRUE(k_0_head <= k.dim(0) && k.dim(0) < k_0_tail);
	EXPECT_TRUE(k_1_head <= k.dim(1) && k.dim(1) < k_1_tail);
      }
      return 0;
    }
  };

  TEST_F(DataStoreViewTest, Map) {
    MapTestMapper0 mapper0(nprocs);
    ds_->map(mapper0, *view_);
  }

#ifdef BACKEND_KMR
  TEST_F(DataStoreViewTest, Collate) {
    // At first, data elements are located on rank 0 only.
    {
      kmrnext::Key key(2);
      for (size_t i = 0; i < ds_->dim(0); i++) {
	key.set_dim(0, i);
	for (size_t j = 0; j < ds_->dim(1); j++) {
	  key.set_dim(1, j);
	  ds_->get(key);
	  EXPECT_EQ(0, ds_->data_element_at(key)->owner());
	}
      }
    }

    ds_->collate();
    if (nprocs > 1) {
      EXPECT_TRUE(ds_->collated());
    }

    // As the default Split is <All, None>, each process has data elements
    // in two rows as a result of data relocation
    {
      int owner;
      kmrnext::Key key(2);
      for (size_t i = 0; i < ds_->dim(0); i++) {
	owner = static_cast<int>(i / 2);
	key.set_dim(0, i);
	for (size_t j = 0; j < ds_->dim(1); j++) {
	  key.set_dim(1, j);
	  ds_->get(key);
	  EXPECT_EQ(owner, ds_->data_element_at(key)->owner());
	}
      }
    }

    ds_->set_split(*split_);
    ds_->collate();
    if (nprocs > 1) {
      EXPECT_TRUE(ds_->collated());
    }

    // As a detail split is set, data elements are grouped and stored on
    // different processes
    {
      int *owners = create_owners();
      kmrnext::Key key(2);
      for (size_t i = 0; i < ds_->dim(0); i++) {
	key.set_dim(0, i);
	for (size_t j = 0; j < ds_->dim(1); j++) {
	  key.set_dim(1, j);
	  ds_->get(key);
	  EXPECT_EQ(owners[j * ds_->dim(0) + i],
	   	    ds_->data_element_at(key)->owner());
	}
      }
      delete_owners(owners);
    }

    ds_->collate();
    if (nprocs > 1) {
      EXPECT_FALSE(ds_->collated());
    }
  }

  class MapTestMapper1 : public kmrnext::DataStore::Mapper {
    int nprocs_;
    int nprocs_task_;
  public:
    MapTestMapper1(int nprocs, int nprocs_task) :
      nprocs_(nprocs), nprocs_task_(nprocs_task) {}

    int operator()(kmrnext::DataStore *inds, kmrnext::DataStore *outds,
		   kmrnext::Key& key, std::vector<kmrnext::DataPack>& dps,
		   kmrnext::DataStore::MapEnvironment& env)
    {
      size_t blk_len0, blk_len1;
      {
	if (nprocs_ > 2 && nprocs_ % 2 == 0) {
	  blk_len0 = inds->dim(0) / env.view.dim(0);
	} else {
	  blk_len0 = inds->dim(0);
	}
	blk_len1 = inds->dim(1) / env.view.dim(1);
      }

      int nprocs_sub;
      MPI_Comm_size(env.mpi_comm, &nprocs_sub);
      EXPECT_EQ(nprocs_task_, nprocs_sub);

      int local_count = static_cast<int>(dps.size());
      int total_count;
      MPI_Allreduce(&local_count, &total_count, 1, MPI_INT, MPI_SUM,
		    env.mpi_comm);
      EXPECT_EQ(blk_len0 * blk_len1, total_count);

      size_t k_0_head, k_0_tail, k_1_head, k_1_tail;
      {
	kmrnext::DataPack dp0 = dps.at(0);
	kmrnext::Key k = dp0.key();
	k_0_head = (k.dim(0) / blk_len0) * blk_len0;
	k_0_tail = k_0_head + blk_len0;
	k_1_head = (k.dim(1) / blk_len1) * blk_len1;
	k_1_tail = k_1_head + blk_len1;
      }
      for (std::vector<kmrnext::DataPack>::iterator itr = dps.begin();
	   itr != dps.end(); itr++) {
	kmrnext::Key k = (*itr).key();
	EXPECT_TRUE(k_0_head <= k.dim(0) && k.dim(0) < k_0_tail);
	EXPECT_TRUE(k_1_head <= k.dim(1) && k.dim(1) < k_1_tail);
      }
      return 0;
    }
  };

  TEST_F(DataStoreViewTest, MapKMR) {
    // As the default Split is <ALL, None, ...>, data relocation occurs.
    int nprocs_task;
    if (nprocs > 2 && nprocs % 2 == 0) {
      nprocs_task = 2;
    } else {
      nprocs_task = nprocs;
    }
    MapTestMapper1 mapper0(nprocs, nprocs_task);
    ds_->map(mapper0, *view_);
    if (nprocs > 1) {
      EXPECT_TRUE(ds_->collated());
    }

    // As the View and Split is same, each task runs on a process.
    nprocs_task = 1;
    MapTestMapper1 mapper1(nprocs, nprocs_task);
    ds_->set_split(*split_);
    ds_->map(mapper1, *view_);
    if (nprocs > 1) {
      EXPECT_TRUE(ds_->collated());
    }

    // As using the same view, data relocation does not happen.
    ds_->map(mapper1, *view_);
    if (nprocs > 1) {
      EXPECT_FALSE(ds_->collated());
    }
  }
#endif

}
