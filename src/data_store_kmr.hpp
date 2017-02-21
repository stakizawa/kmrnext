#ifndef DATA_STORE_KMR_HPP
#define DATA_STORE_KMR_HPP
/// \file
/// DataStore implementation header for KMR backend

/// The current implementation of KMR-based DataStore consumes huge amount
/// of memory.  It uses just an array, dense matrix, to store Data, but
/// actually it is sparse.  It should be modified.
///
/// The basic storategy for parallelization
/// add        : Rank 0 process adds data locally in case of serial context.
///              Each process adds data locally in case of parallel context
///              (when using map() function).
/// get        : A process that has the data broadcasts it.
/// get<view>  : The data is allgathered between all processes.
/// map        : Before applying the fanctor, the data elements are shuffled
///              between nodes according to the specified split pattern.
///              When applying the fanctor, the DataStore is split by the
///              specified (task) view, which defines a unit of task, nodes
///              that have the data elements which should be processed as a
///              task are grouped, and then the functor runs on the nodes
///              to process the task.
/// load_xxx   : Array elements are added on rank 0 only, but they are
///              scattered when the data in them are loaded.  So, if there
///              are enough number of nodes (#nodes > #array elements), the
///              array elements can be simultaneously loaded.

#include "data_store.hpp"
#include "data_element.hpp"
#include <algorithm>

// Count of vector elements pre-reserved for adding elements in OpenMP
// critical region.
#define CRITICAL_VECTOR_PRE_RESERVE_SIZE 100

// If 1, it performs collate() using KMR
#define COLLATE_KMR 0

namespace kmrnext{
  using namespace std;

  /////////////////////////////////////////////////////////////////////////
  // Basic implementation of DataStore
  /////////////////////////////////////////////////////////////////////////
  template <typename DE>
  class DataStoreImpl : public DataStoreCommonImpl {
    typedef DataStoreCommonImpl base;

  public:

    virtual ~DataStoreImpl();

    void set(const size_t* val);

    void set_dim(const size_t idx, const size_t siz);

    virtual void add(const Key& key, const Data& data);

    virtual DataPack get(const Key& key);

    virtual vector<DataPack>* get(const View& view, const Key& key);

    virtual DataPack remove(const Key& key);

    // The implementation is common to both Serial and KMR backend
    virtual void set_from(const vector<DataStore*>& dslist);

    // The implementation is common to both Serial and KMR backend
    virtual void split_to(vector<DataStore*>& dslist);

    virtual void map(Mapper& m, const View& view, DataStore* outds=self_);

    void set_split(const View& split);

    View get_split();

    virtual void collate();

    bool collated();

    void set_force_collate(bool status);

    // The implementation is common to both Serial and KMR backend
    void load_files(const vector<string>& files, Loader<string>& loader);

    // The implementation is common to both Serial and KMR backend
    void load_integers(const vector<long>& ints, Loader<long>& loader);

    void load_parallel(Loader<long>& loader);

    // The implementation is common to both Serial and KMR backend
    virtual DataStore* duplicate();

    string dump(DataPack::Dumper& dumper);

    long count();

    size_t key_to_split_index(const Key& key, const View& view);

    // The implementation is common to both Serial and KMR backend
    DataElement& data_element_at(const Key& key);

  protected:
    // Stored DataElements
    vector<DE> dlist_;
    // Size of dlist_
    size_t dlist_size_;
    // True if the input and output DataStore of map function is same
    bool map_inplace_;
    // A KMRNext object that stores execution status
    KMRNext* kmrnext_;
    // True if the DataStore should be processed in parallel
    bool parallel_;
    // Split pattern of the DataStore, that defines data distribution
    View* split_;
    // Set to be true if the last call of map() or collate() actually
    // performed collate.
    bool collated_;
    // If true, collating data is always performed when map() or collate()
    // is called.
    bool force_collate_;

    // It creates a new DataStore.
    //
    // Size of each dimension should be set by DataStore.set() function
    // before using this DataStore.
    //
    // \param[in] size number of dimensions of this DataStore.
    // \param[in] kn   an instance of KMRNext context
    DataStoreImpl(const size_t siz, KMRNext* kn);

  };

  /////////////////////////////////////////////////////////////////////////
  // An implementation of DataStore that stores data in memory
  /////////////////////////////////////////////////////////////////////////
  class InMemoryDataStore : public DataStoreImpl<DataElement> {
    typedef DataStoreImpl<DataElement> base;

  public:
    /// It creates a new InMemoryDataStore.
    ///
    /// Size of each dimension should be set by DataStore.set() function
    /// before using this DataStore.
    ///
    /// \param[in] size number of dimensions of this DataStore.
    /// \param[in] kn   an instance of KMRNext context
    InMemoryDataStore(const size_t siz, KMRNext* kn) : base(siz, kn) {}

    virtual ~InMemoryDataStore() {}
  };

  /////////////////////////////////////////////////////////////////////////
  // An implementation of DataStore that stores data in files on processes
  /////////////////////////////////////////////////////////////////////////
  class SimpleFileDataStore : public DataStoreImpl<SimpleFileDataElement> {
    typedef DataStoreImpl<SimpleFileDataElement> base;

  public:
    /// It creates a new SimpleFileDataStore.
    ///
    /// Size of each dimension should be set by DataStore.set() function
    /// before using this DataStore.
    ///
    /// \param[in] size number of dimensions of this DataStore.
    /// \param[in] kn   an instance of KMRNext context
    SimpleFileDataStore(const size_t siz, KMRNext* kn)
      : base(siz, kn), data_updated_(false), data_cached_(false) {}

    // The implementation is common to both Serial and KMR backend
    virtual ~SimpleFileDataStore();

    // The implementation is common to both Serial and KMR backend
    void add(const Key& key, const Data& data);

    DataPack get(const Key& key);

    vector<DataPack>* get(const View& view, const Key& key);

    DataPack remove(const Key& key);

    // The implementation is common to both Serial and KMR backend
    void set_from(const vector<DataStore*>& dslist);

    // The implementation is common to both Serial and KMR backend
    void split_to(vector<DataStore*>& dslist);

    // The implementation is common to both Serial and KMR backend
    void map(Mapper& m, const View& view, DataStore* outds=self_);

    // The implementation is common to both Serial and KMR backend
    DataStore* duplicate();

    virtual void collate();

  private:
    // True if the data in DataStore is updated, but not written to a file
    bool data_updated_;
    // True if the Data in a file is cached in memory
    bool data_cached_;

    // It returns name of file that stores data elements of the DataStore.
    string filename();

    // It writes data elements in the DataStore to a file.
    // The implementation is common to both Serial and KMR backend
    //
    // \return   true if Data are stored in a file
    bool store();

    // It reads data elements in a file and then loads to the DataStore.
    // The implementation is common to both Serial and KMR backend
    //
    // \return   true if Data are loaded from a file
    bool load();

    // It clears memory cache.
    // The implementation is common to both Serial and KMR backend
    void clear_cache();

  };

}

//##########################################################################//
// Method implementations specialized for the KMR backend
//##########################################################################//

//////////////////////////////////////////////////////////////////////////////
// This namespace contains functions and structures used to implement
// KMR backend functions
//////////////////////////////////////////////////////////////////////////////
namespace kmrutils {
  using namespace std;
  using namespace kmrnext;

  ////////////////////////////////////////////////////////////////////////////
  // Parameter class and task for mapper_get_view
  //
  // The task copies a key-value to vector<DataPack>.
  //
  // The task is a KMR mapper function.  Key of the key-value should be an
  // integer and value of the key-value should be a DataPack.  It also
  // sets shared flags to data in the DataStore.
  ////////////////////////////////////////////////////////////////////////////
  template <typename DE>
  struct param_mapper_get_view {
    DataStoreImpl<DE>* ds;    // target DataStore
    vector<DE>&        dlist; // DataElements in the DataStore
    vector<DataPack>*  dps;   // result vector
  };
  template <typename DE>
  int mapper_get_view(const struct kmr_kv_box kv0, const KMR_KVS* kvi,
		      KMR_KVS* kvo, void* p, const long i) {
    param_mapper_get_view<DE>* param = (param_mapper_get_view<DE>*)p;
    size_t idx = kv0.k.i;
    Key key = param->ds->index_to_key(idx);
    if (!param->dlist[idx].is_set()) {
      int owner;
      memcpy(&owner, kv0.v.p, sizeof(int));
      param->dlist[idx].set(const_cast<char*>(kv0.v.p) + sizeof(int),
			    kv0.vlen - sizeof(int));
      param->dlist[idx].set_owner(owner);
    }
    param->dlist[idx].shared();
    Data dat(param->dlist[idx].value(), param->dlist[idx].size());
#if KMR_OMP
    // As a KMR map function can be run in parallel by OpenMP, shared resources
    // should be modified in critical regions.
#ifdef _OPENMP
    #pragma omp critical
#endif
#endif
    param->dps->push_back(DataPack(key, dat, true));
    return MPI_SUCCESS;
  }

  ////////////////////////////////////////////////////////////////////////////
  // Parameter class and task for mapper_map
  //
  // Tha task maps Data in the input DataStore.
  //
  // Tha task is a KMR mapper function.  Key and value of the key-value
  // should be an integer that represents index of Key-Data vactor.
  ////////////////////////////////////////////////////////////////////////////
  template <typename DE>
  struct param_mapper_map {
    DataStore::Mapper&          mapper;
    DataStoreImpl<DE>*          ids;
    DataStore*                  ods;
    const View&                 view;
    DataStore::MapEnvironment&  env;
    vector< vector<DataPack> >& dpgroups;  // Data to be processed
    bool                        map_local;
  };
  template <typename DE>
  int mapper_map(const struct kmr_kv_box kv0, const KMR_KVS* kvi,
		 KMR_KVS* kvo, void* p, const long i) {
    param_mapper_map<DE>* param = static_cast<param_mapper_map<DE>*>(p);
    size_t idx = kv0.k.i;
    vector<DataPack>& dps = param->dpgroups.at(idx);
    if (dps.size() != 0) {
      Key viewed_key =
	param->ids->key_to_viewed_key(dps.at(0).key(), param->view);
      if (!param->map_local) {
	param->env.mpi_comm = kvi->c.mr->comm;
      }
      param->mapper(param->ids, param->ods, viewed_key, dps, param->env);
    }
    return MPI_SUCCESS;
  }

  // It is a wrapper of KMR map function that locally maps key-values in
  // the input KVS.
  template <typename DE>
  int map_local(KMR_KVS* kvi, KMR_KVS* kvo, void* arg, kmr_mapfn_t m) {
    param_mapper_map<DE>* param = static_cast<param_mapper_map<DE>*>(arg);
    MPI_Comm self_comm;
    MPI_Comm_split(kvi->c.mr->comm, kvi->c.mr->rank, kvi->c.mr->rank,
		   &self_comm);
    param->env.mpi_comm = self_comm;
    struct kmr_option kmr_nothrd = kmr_noopt;
    kmr_nothrd.nothreading = 1;
    kmr_map(kvi, kvo, arg, kmr_nothrd, m);
    MPI_Comm_free(&self_comm);
    return MPI_SUCCESS;
  }

  // A data set used in collate()
  class CollatePack {
  public:
    Key key_;
    DataElement* de_;
    CollatePack(Key& k, DataElement* de) : key_(k), de_(de) {}
  };

  // It performs collate() using KMR function.
  void collate_kmr(KMRNext *kn, DataStore* ds,
		   vector< vector<CollatePack> >& cpgroups,
		   size_t ncpgroups);

  // It performs collate() using MPI alltoall function.
  void collate_mpi(KMRNext *kn, DataStore* ds,
		   vector< vector<CollatePack> >& cpgroups,
		   size_t ncpgroups);

}

//////////////////////////////////////////////////////////////////////////////
// This namespace contains class implementations for the KMR backend
//////////////////////////////////////////////////////////////////////////////
namespace kmrnext {
  using namespace std;
  using namespace kmrutils;

  template <typename DE>
  DataStoreImpl<DE>::DataStoreImpl(const size_t siz, KMRNext* kn)
    : base(siz), dlist_(vector<DE>()), dlist_size_(0), map_inplace_(false),
      kmrnext_(kn), parallel_(false), split_(NULL), collated_(false),
      force_collate_(false) {}

  template <typename DE>
  DataStoreImpl<DE>::~DataStoreImpl() {
    if (split_ != NULL) {
      delete split_;
    }
  }

  template <typename DE>
  void DataStoreImpl<DE>::add(const Key& key, const Data& data) {
#if VALIDATION
    check_key_range(key);
#endif
    if (dlist_size_ == 0) {
      set(value_);
    }
    if (parallel_ || kmrnext_->rank() == 0) {
      size_t idx = key_to_index(key);
      try {
	if (map_inplace_) {
	  dlist_[idx].replace(data.value(), data.size());
	} else {
	  dlist_[idx].set(data.value(), data.size());
	}
      }
      catch (runtime_error& e) {
	cerr << "Failed to add a data to DataStore" << to_string()
	     << " at Key" << key.to_string()
	     << "on Rank" << kmrnext_->rank() << "." << endl;
	throw e;
      }
      dlist_[idx].set_owner(kmrnext_->rank());
    }
  }

  template <typename DE>
  DataPack DataStoreImpl<DE>::get(const Key& key) {
#if VALIDATION
    check_key_range(key);
#endif

    size_t idx = key_to_index(key);
    if (dlist_[idx].is_shared()) {
      Data dat(dlist_[idx].value(), dlist_[idx].size());
      return DataPack(key, dat, true);
    }

    KMR_KVS* snd = kmr_create_kvs(kmrnext_->kmr(),
				  KMR_KV_INTEGER, KMR_KV_OPAQUE);
    KMR_KVS* rcv = kmr_create_kvs(kmrnext_->kmr(),
				  KMR_KV_INTEGER, KMR_KV_OPAQUE);
    if (dlist_[idx].is_set()) {
      size_t snd_siz = sizeof(int) + dlist_[idx].size();
      char* snd_buf = new char[snd_siz];
      int owner = kmrnext_->rank();
      memcpy(snd_buf, &owner, sizeof(int));
      memcpy(snd_buf + sizeof(int), dlist_[idx].value(), dlist_[idx].size());
      struct kmr_kv_box kv;
      kv.klen = (int)sizeof(long);
      kv.vlen = (int)snd_siz;
      kv.k.i  = idx;
      kv.v.p  = snd_buf;
      kmr_add_kv(snd, kv);
      delete[] snd_buf;
    }
    kmr_add_kv_done(snd);
    kmr_replicate(snd, rcv, kmr_noopt);
    long local_count;
    kmr_local_element_count(rcv, &local_count);
#ifdef DEBUG
    if (local_count > 1) {
      kmr_free_kvs(rcv);
      throw runtime_error("There are too many data.");
    }
#endif
    if (local_count == 1) {
      struct kmr_kv_box kv;
      kmr_take_one(rcv, &kv);
      if (!dlist_[idx].is_set()) {
	// Copy data for future access.
	int owner;
	memcpy(&owner, kv.v.p, sizeof(int));
	dlist_[idx].set((const_cast<char*>(kv.v.p) + sizeof(int)),
			kv.vlen - sizeof(int));
	dlist_[idx].set_owner(owner);
      }
      dlist_[idx].shared();
    }
    kmr_free_kvs(rcv);

    Data dat(dlist_[idx].value(), dlist_[idx].size());
    return DataPack(key, dat, true);
  }

  template <typename DE>
  vector<DataPack>* DataStoreImpl<DE>::get(const View& view, const Key& key) {
#if VALIDATION
    check_view(view);
    check_key_range(key);
#endif
    vector<DataPack>* dps = new vector<DataPack>();

    size_t* blk_sizs = new size_t[size_];
    for (size_t i = 0; i < size_; i++) {
      blk_sizs[i] = (view.dim(i) > 0)? (value_[i] / view.dim(i)) : 1;
    }

    KMR_KVS* snd = kmr_create_kvs(kmrnext_->kmr(),
				  KMR_KV_INTEGER, KMR_KV_OPAQUE);
    KMR_KVS* rcv = kmr_create_kvs(kmrnext_->kmr(),
				  KMR_KV_INTEGER, KMR_KV_OPAQUE);
#if KMR_OMP
#ifdef _OPENMP
    #pragma omp parallel for schedule(static, OMP_FOR_CHUNK_SIZE)
#endif
#endif
    for (size_t i = 0; i < dlist_size_; i++) {
      if (!dlist_[i].is_set()) {
	continue;
      }
      Key tmpkey = index_to_key(i);
      bool match = true;
      for (size_t j = 0; j < size_; j++) {
	if ((view.dim(j) == View::SplitAll && key.dim(j) != tmpkey.dim(j)) ||
	    (view.dim(j) > 0 && key.dim(j) != tmpkey.dim(j) / blk_sizs[j])) {
	  match = false;
	  break;
	}
      }
      if (match) {
	if (dlist_[i].is_shared()) {
	  // the data is already replicated
	  Data dat(dlist_[i].value(), dlist_[i].size());
#if KMR_OMP
#ifdef _OPENMP
          #pragma omp critical
#endif
#endif
	  dps->push_back(DataPack(tmpkey, dat, true));
	  continue;
	} else {
	  // replicate the data
	  size_t snd_siz = sizeof(int) + dlist_[i].size();
	  char* snd_buf = new char[snd_siz];
	  int owner = kmrnext_->rank();
	  memcpy(snd_buf, &owner, sizeof(int));
	  memcpy(snd_buf + sizeof(int), dlist_[i].value(), dlist_[i].size());
	  struct kmr_kv_box kv;
	  kv.klen = static_cast<int>(sizeof(size_t));
	  kv.vlen = static_cast<int>(snd_siz);
	  kv.k.i  = i;
	  kv.v.p  = snd_buf;
	  kmr_add_kv(snd, kv);
	  delete[] snd_buf;
	}
      }
    }
    kmr_add_kv_done(snd);
    kmr_replicate(snd, rcv, kmr_noopt);
    param_mapper_get_view<DE> param = { this, dlist_, dps };
    kmr_map(rcv, NULL, (void*)&param, kmr_noopt, mapper_get_view<DE>);

    delete[] blk_sizs;
    return dps;
  }

  template <typename DE>
  DataPack DataStoreImpl<DE>::remove(const Key& key) {
    DataPack dp0 = get(key);
    if (dp0.data().value() != NULL) {
      size_t idx = key_to_index(key);
      dlist_[idx].clear();
    }
    return dp0;
  }

  template <typename DE>
  void DataStoreImpl<DE>::map(Mapper& m, const View& view, DataStore* outds) {
#if VALIDATION
    check_map_args(view, outds);
#endif
    if (dlist_size_ == 0) {
      return;
    }

    collate();
    View split = get_split();

    size_t nkeys = 1;
    bool is_local = true;
    {
      bool is_same_view = true;
      for (size_t i = 0; i < size_; i++) {
	if (view.dim(i) == View::SplitAll) {
	  nkeys *= value_[i];
	} else if (view.dim(i) > 0) { // split count is specified
	  nkeys *= view.dim(i);
	}
	if (view.dim(i) != split.dim(i)) {
	  is_same_view = false;
	}
      }
      if (!is_same_view) {
	is_local = false;
      }
    }

    vector< vector<DataPack> > dpgroups(nkeys);
#ifdef _OPENMP
    #pragma omp parallel
    {
      vector< vector<DataPack> > dpgroups_p(nkeys);
      #pragma omp for schedule(static, OMP_FOR_CHUNK_SIZE)
      for (size_t i = 0; i < dlist_size_; i++) {
	if (!dlist_[i].is_set() ||
	    (dlist_[i].is_shared() && dlist_[i].owner() != kmrnext_->rank())) {
	  continue;
	}
	Key tmpkey = index_to_key(i);
	size_t viewed_idx = key_to_viewed_index(tmpkey, view);
	vector<DataPack>& dps = dpgroups_p.at(viewed_idx);
	Data dat(dlist_[i].value(), dlist_[i].size());
	dps.push_back(DataPack(tmpkey, dat));
      }
      #pragma omp critical
      {
	for (size_t i = 0; i < nkeys; i++) {
	  vector<DataPack>& dps_p = dpgroups_p.at(i);
	  if (dps_p.size() > 0) {
	    vector<DataPack>& dps = dpgroups.at(i);
	    if (dps.size() == 0) {
	      dps.reserve(CRITICAL_VECTOR_PRE_RESERVE_SIZE);
	    }
	    copy(dps_p.begin(), dps_p.end(), back_inserter(dps));
	  }
	}
      }
    }
#else
    for (size_t i = 0; i < dlist_size_; i++) {
      if (!dlist_[i].is_set() ||
	  (dlist_[i].is_shared() && dlist_[i].owner() != kmrnext_->rank())) {
	continue;
      }
      Key tmpkey = index_to_key(i);
      size_t viewed_idx = key_to_viewed_index(tmpkey, view);
      vector<DataPack>& dps = dpgroups.at(viewed_idx);
      Data dat(dlist_[i].value(), dlist_[i].size());
      dps.push_back(DataPack(tmpkey, dat));
    }
#endif

    if (kmrnext_->profile()) {
      long dlist_count = 0;
      for (size_t i = 0; i < nkeys; i++) {
	vector<DataPack> &dps = dpgroups.at(i);
	if (dps.size() > 0) {
	  dlist_count += 1;
	}
      }
      ostringstream os;
      os << "count of data to be mapped: " << dlist_count;
      profile_out(kmrnext_, os.str());
    }

    KMR_KVS* ikvs = kmr_create_kvs(kmrnext_->kmr(),
				   KMR_KV_INTEGER, KMR_KV_INTEGER);
#if KMR_OMP
#ifdef _OPENMP
    #pragma omp parallel for
#endif
#endif
    for (size_t i = 0; i < nkeys; i++) {
      vector<DataPack> &dps = dpgroups.at(i);
      if (dps.size() > 0) {
	struct kmr_kv_box kv;
	kv.klen = static_cast<int>(sizeof(size_t));
	kv.vlen = static_cast<int>(sizeof(size_t));
	kv.k.i  = i;
	kv.v.i  = i;
	kmr_add_kv(ikvs, kv);
      }
    }
    kmr_add_kv_done(ikvs);

    DataStoreImpl<DE>* _outds = static_cast<DataStoreImpl<DE>*>(outds);
    if (outds == self_ || outds == this) {
      map_inplace_ = true;
      _outds = this;
    }
    _outds->parallel_ = true;
    MapEnvironment env = { kmrnext_->rank(), view, split,
			   MPI_COMM_NULL };
    param_mapper_map<DE> param = { m, static_cast<DataStoreImpl<DE>*>(this),
				   _outds, view, env, dpgroups, is_local };
    if (is_local) {
      double timestamp[2] = {0, 0};
      long task_count;
      kmr_local_element_count(ikvs, &task_count);
      if (kmrnext_->profile()) {
	timestamp[0] = gettime(kmrnext_->kmr()->comm);
      }

      map_local<DE>(ikvs, NULL, (void*)&param, mapper_map<DE>);

      if (kmrnext_->profile()) {
	timestamp[1] = gettime(kmrnext_->kmr()->comm);
	double total_time = (timestamp[1] - timestamp[0]) / 1E9;
	ostringstream os;
	os << "timing of map locally: total_exec="
	   << fixed << total_time << " each_exec="
	   << fixed << (total_time / static_cast<double>(task_count))
	   << "(sec), "
	   << "count of task execution: " << task_count;
	profile_out(kmrnext_, os.str());
      }
    } else {
      kmr_map_multiprocess_by_key(ikvs, NULL, static_cast<void*>(&param),
				  kmr_noopt, kmrnext_->rank(), mapper_map<DE>);
    }
    _outds->parallel_ = false;
    if (outds == self_ || outds == this) {
      map_inplace_ = false;
      // unshare all data
#ifdef _OPENMP
      #pragma omp parallel for schedule(static, OMP_FOR_CHUNK_SIZE)
#endif
      for (size_t i = 0; i < dlist_size_; i++) {
	if (dlist_[i].is_shared() && dlist_[i].owner() != kmrnext_->rank()) {
	  dlist_[i].clear();
	  continue;
	}
	dlist_[i].unshared();
      }
    }
  }

  template <typename DE>
  void DataStoreImpl<DE>::set_split(const View& split) {
#if VALIDATION
    check_view(split);
#endif
    if (split_ == NULL) {
      split_ = new View(split.size());
    }
    for (size_t i = 0; i < split.size(); i++) {
      split_->set_dim(i, split.dim(i));
    }
  }

  template <typename DE>
  View DataStoreImpl<DE>::get_split() {
    if (split_ == NULL) {
      split_ = new View(size_);
      for (size_t i = 0; i < size_; i++) {
	if (i == 0) {
	  split_->set_dim(i, View::SplitAll);
	} else {
	  split_->set_dim(i, View::SplitNone);
	}
      }
    }
    return *split_;
  }

  template <typename DE>
  void DataStoreImpl<DE>::collate() {
    double time_collate_start = 0, time_collate_finish = 0;
    if (kmrnext_->profile()) {
      time_collate_start = gettime(kmrnext_->kmr()->comm);
    }

    // Split pattern set to this DataStore
    View split = get_split();
    // Count of viewed data
    size_t ndata;
    // Count of viewed data this process should have
    size_t range_cnt;
    // The first index of viewed data this process should have
    size_t range_start = 0;
    // The last index of viewed data this process should have
    size_t range_end   = 0;
    // Count of elements in each split data block
    size_t each_cnt = 1;

    // Calculate parameters
    {
      ndata = 1;
      for (size_t i = 0; i < size_; i++) {
	if (split.dim(i) == View::SplitAll) {
	  ndata *= value_[i];
	} else if (split.dim(i) == View::SplitNone) {
	  each_cnt *= value_[i];
	} else {
	  // Split count is specified
#ifdef DEBUG
	  assert(value_[i] % split.dim(i) == 0);
#endif
	  ndata *= split.dim(i);
	  each_cnt *= (value_[i] / split.dim(i));
	}
      }

      // Calculate indices of Data this process should hold
      {
	size_t avg_cnt = ndata / static_cast<size_t>(kmrnext_->nprocs());
	size_t rem_cnt = ndata % static_cast<size_t>(kmrnext_->nprocs());
	range_cnt =
	  (kmrnext_->rank() < static_cast<int>(rem_cnt))?
	  avg_cnt + 1 : avg_cnt;
	if (range_cnt != 0) {
	  range_start = avg_cnt * static_cast<size_t>(kmrnext_->rank())
	    + ((kmrnext_->rank() < static_cast<int>(rem_cnt))?
	       kmrnext_->rank() : rem_cnt);
	  range_end = range_start + range_cnt - 1;
	}
      }
    }

    // Check if collate is necessary.
    if (force_collate_) {
      collated_ = true;
    } else {
      // It holds the DataPacks this process should hold
      vector<int> data_cnts(range_cnt, 0);
      if (range_cnt != 0) {
#ifdef _OPENMP
        #pragma omp parallel
	{
	  vector<int> data_cnts_p(range_cnt, 0);
          #pragma omp for schedule(static, OMP_FOR_CHUNK_SIZE)
	  for (size_t i = 0; i < dlist_size_; i++) {
	    if (!dlist_[i].is_set() ||
		(dlist_[i].is_shared() &&
		 dlist_[i].owner() != kmrnext_->rank())) {
	      continue;
	    }
	    Key tmpkey = index_to_key(i);
	    size_t viewed_idx = key_to_split_index(tmpkey, split);
	    if (range_start <= viewed_idx && viewed_idx <= range_end) {
	      size_t idx = viewed_idx - range_start;
	      data_cnts_p[idx] = data_cnts_p[idx] + 1;
	    }
	  }
          #pragma omp critical
	  {
	    for (size_t i = 0; i < range_cnt; i++) {
	      data_cnts[i] = data_cnts[i] + data_cnts_p[i];
	    }
	  }
	}
#else
	for (size_t i = 0; i < dlist_size_; i++) {
	  if (!dlist_[i].is_set() ||
	      (dlist_[i].is_shared() &&
	       dlist_[i].owner() != kmrnext_->rank())) {
	    continue;
	  }
	  Key tmpkey = index_to_key(i);
	  size_t viewed_idx = key_to_split_index(tmpkey, split);
	  if (range_start <= viewed_idx && viewed_idx <= range_end) {
	    size_t idx = viewed_idx - range_start;
	    data_cnts[idx] = data_cnts[idx] + 1;
	  }
	}
#endif
      }

      int local_need_collate = 0;
      for (size_t i = 0; i < range_cnt; i++) {
	if (data_cnts[i] != each_cnt) {
	  local_need_collate = 1;
	  break;
	}
      }
      int need_collate;
      MPI_Allreduce(&local_need_collate, &need_collate, 1, MPI_INT, MPI_MAX,
		    kmrnext_->kmr()->comm);
      if (kmrnext_->profile()) {
	if (kmrnext_->rank() == 0) {
	  ostringstream os;
	  if (need_collate == 0) {
	    os << "No need to collate.";
	  } else {
	    os << "Perform collate.";
	  }
	  profile_out(kmrnext_, os.str());
	}
      }

      if (need_collate == 0) {
	// no need to collate
	collated_ = false;
	return;
      } else {
	collated_ = true;
      }
    }

    // Collate from here

    // Search data should be sent
    vector< vector<CollatePack> > cpgroups(ndata);
    {
#ifdef _OPENMP
      #pragma omp parallel
      {
	vector< vector<CollatePack> > cpgroups_p(ndata);
        #pragma omp for schedule(static, OMP_FOR_CHUNK_SIZE)
	for (size_t i = 0; i < dlist_size_; i++) {
	  if (!dlist_[i].is_set()) {
	    continue;
	  }
	  if (dlist_[i].is_shared() && dlist_[i].owner() != kmrnext_->rank()) {
	    dlist_[i].clear();
	    continue;
	  }
	  dlist_[i].unshared();
	  Key tmpkey = index_to_key(i);
	  size_t viewed_idx = key_to_split_index(tmpkey, split);
	  if (range_cnt > 0 &&
	      !(range_start <= viewed_idx && viewed_idx <= range_end)) {
	    vector<CollatePack>& cps = cpgroups_p.at(viewed_idx);
	    cps.push_back(CollatePack(tmpkey, &dlist_[i]));
	  }
	}
        #pragma omp critical
	{
	  for (size_t i = 0; i < ndata; i++) {
	    vector<CollatePack>& cps_p = cpgroups_p.at(i);
	    if (cps_p.size() > 0) {
	      vector<CollatePack>& cps = cpgroups.at(i);
	      if (cps.size() == 0) {
		cps.reserve(CRITICAL_VECTOR_PRE_RESERVE_SIZE);
	      }
	      copy(cps_p.begin(), cps_p.end(), back_inserter(cps));
	    }
	  }
	}
      }
#else
      for (size_t i = 0; i < dlist_size_; i++) {
	if (!dlist_[i].is_set()) {
	  continue;
	}
	if (dlist_[i].is_shared() && dlist_[i].owner() != kmrnext_->rank()) {
	  dlist_[i].clear();
	  continue;
	}
	dlist_[i].unshared();
	Key tmpkey = index_to_key(i);
	size_t viewed_idx = key_to_split_index(tmpkey, split);
	if (range_cnt > 0 &&
	    !(range_start <= viewed_idx && viewed_idx <= range_end)) {
	  vector<CollatePack>& cps = cpgroups.at(viewed_idx);
	  cps.push_back(CollatePack(tmpkey, &dlist_[i]));
	}
      }
#endif
    }

#if COLLATE_KMR
    collate_kmr(kmrnext_, this, cpgroups, ndata);
#else
    collate_mpi(kmrnext_, this, cpgroups, ndata);
#endif

    if (kmrnext_->profile()) {
      time_collate_finish = gettime(kmrnext_->kmr()->comm);
      if (kmrnext_->rank() == 0) {
	ostringstream os;
	os << "Time spent for collate: "
	   << (time_collate_finish - time_collate_start) << " nano-sec";
	profile_out(kmrnext_, os.str());
      }
    }
  }

  template <typename DE>
  bool DataStoreImpl<DE>::collated() {
    return collated_;
  }

  template <typename DE>
  void DataStoreImpl<DE>::set_force_collate(bool status) {
    force_collate_ = status;
  }

  template <typename DE>
  void DataStoreImpl<DE>::load_parallel(Loader<long>& loader) {
#if 1
    // load in parallel without using map().

    parallel_ = true;
    loader(this, kmrnext_->rank());
    parallel_ = false;
#else
    // load in parallel using map().

    // Create a 1D DataStore that stores ranks
    DataStoreImpl<DE> ds0(1, kmrnext_);
    ds0.set_dim(0, kmrnext_->nprocs());
    Key key(1);
    for (long i = 0; i < kmrnext_->nprocs(); i++) {
      key.set_dim(0, i);
      Data dat(&i, sizeof(long));
      ds0.add(key, dat);
    }

    // Use Split <T> so that each process loads an array element.
    View split_ds0(1);
    split_ds0.set_dim(0, View::SplitAll);
    ds0.set_split(split_ds0);

    // Define a mapper for the loader
    class WrappedLoader : public DataStore::Mapper {
    public:
      DataStore::Loader<long>& loader_;
      KMRNext* kmrnext_;

      WrappedLoader(DataStore::Loader<long>& ldr, KMRNext* next)
	: loader_(ldr), kmrnext_(next) {}
      int operator()(DataStore* inds, DataStore* outds,
		     Key& k, vector<DataPack>& dps,
		     DataStore::MapEnvironment& env)
      {
#ifdef DEBUG
	if (dps.size() != 1) {
	  throw runtime_error("System error: Unexpected shuffle behavior.");
	}
#endif
	long* dsval = static_cast<long*>(dps[0].data().value());
#ifdef DEBUG
	if (*dsval != kmrnext_->rank()) {
	  throw runtime_error("System error: Unexpected shuffle behavior.");
	}
#endif
	loader_(outds, *dsval);
	return 0;
      }
    } wloader(loader, kmrnext_);

    // Run the mapper
    View v(1);
    v.set_dim(0, View::SplitAll);
    ds0.map(wloader, v, this);
#endif

    // Set the default split
    View split_ds(size_);
    long nproc_dim = -1;
    {
      for (size_t i = 0; i < size_; i++) {
	split_ds.set_dim(i, View::SplitNone);
      }
      for (size_t i = 0; i < size_; i++) {
	if (value_[i] == static_cast<size_t>(kmrnext_->nprocs())) {
	  nproc_dim = i;
	  break;
	}
      }
    }
    if (nproc_dim != -1) {
      // Use Split <None*, All, None*>, where All is the dimension in the
      // DataStore whose size is equals to the number of processes, so that
      // data elements will be consumed without distribution when the first
      // map() is called.
      split_ds.set_dim(nproc_dim, View::SplitAll);
    } else {
      // Use Split <All, None, ..> and data elements will be split and
      // stored on processes by the first dimension
      split_ds.set_dim(0, View::SplitAll);
    }
    set_split(split_ds);
  }

  template <typename DE>
  string DataStoreImpl<DE>::dump(DataPack::Dumper& dumper) {
    class WrappedDumper : public Mapper {
    public:
      string result_;
      DataPack::Dumper& dumper_;

      WrappedDumper(DataPack::Dumper& dmpr) : dumper_(dmpr) {}
      int operator()(DataStore* inds, DataStore* outds,
		     Key& key, vector<DataPack>& dps,
		     MapEnvironment& env)
      {
	ostringstream os;
	for (vector<DataPack>::iterator itr = dps.begin(); itr != dps.end();
	     itr++) {
	  os << dumper_(*itr);
	}
	string dumped = os.str();
	int local_len = static_cast<int>(dumped.size());
	char* local_cstr = const_cast<char*>(dumped.c_str());
	int nprocs;
	MPI_Comm_size(env.mpi_comm, &nprocs);
	int* local_lens = new int[nprocs];
	MPI_Allgather(&local_len, 1, MPI_INT, local_lens, 1, MPI_INT,
		      env.mpi_comm);
	int total_len = 0;
	int* displs = new int[nprocs];
	for (int i = 0; i < nprocs; i++) {
	  displs[i] = total_len;
	  total_len += local_lens[i];
	}
	total_len += 1; // +1 for '\0'
	char* total_cstr = new char[total_len];
	MPI_Allgatherv(local_cstr, local_len, MPI_CHAR,
		       total_cstr, local_lens, displs, MPI_CHAR, env.mpi_comm);
	total_cstr[total_len - 1] = '\0';

	result_.append(total_cstr);
	delete[] total_cstr;
	delete[] displs;
	delete[] local_lens;
	return 0;
      }
    } dmpr(dumper);

    View view(size_);
    for (size_t i = 0; i < size_; i++) {
      view.set_dim(i, View::SplitNone);
    }
    map(dmpr, view);
    // find master
    int token = (dmpr.result_.size() > 0)? kmrnext_->rank() : -1;
    int master;
    MPI_Allreduce(&token, &master, 1, MPI_INT, MPI_MAX, kmrnext_->kmr()->comm);
    master = (master == -1)? 0 : master;
    // bcast string
    int length = static_cast<int>(dmpr.result_.size()) + 1;
    MPI_Bcast(&length, 1, MPI_INT, master, kmrnext_->kmr()->comm);
    char* result_cstr = new char[length];
    if (kmrnext_->rank() == master) {
      memcpy(result_cstr, dmpr.result_.c_str(), sizeof(char) * length);
    }
    MPI_Bcast(result_cstr, length, MPI_CHAR, master, kmrnext_->kmr()->comm);
    string result(result_cstr);
    delete[] result_cstr;
    return result;
  }

  template <typename DE>
  long DataStoreImpl<DE>::count() {
    class Counter : public Mapper {
    public:
      long result_;
      Counter() : result_(0) {}
      int operator()(DataStore* inds, DataStore* outds,
		     Key& key, vector<DataPack>& dps,
		     MapEnvironment& env)
      {
	long local_count = dps.size();
	MPI_Allreduce(&local_count, &result_, 1, MPI_LONG, MPI_SUM,
	 	      env.mpi_comm);
	return 0;
      }
    } counter;

    View view(size_);
    for (size_t i = 0; i < size_; i++) {
      view.set_dim(i, View::SplitNone);
    }
    map(counter, view);
    long result;
    MPI_Allreduce(&counter.result_, &result, 1, MPI_LONG, MPI_MAX,
		  kmrnext_->kmr()->comm);
    return result;
  }

  template <typename DE>
  size_t DataStoreImpl<DE>::key_to_split_index(const Key& key,const View& view)
  {
    // The implementation is the same as key_to_viewed_index().
#if 1
    // It returns column-ordered index
    size_t idx = 0;
    for (long i = static_cast<long>(size_) - 1; i >= 0; i--) {
      if (view.dim(i) == View::SplitAll || view.dim(i) != View::SplitNone) {
	size_t key_idx = key.dim(i);
	if (view.dim(i) > 0) { // split count is specified
	  size_t blk_i = value_[i] / view.dim(i);
	  key_idx = key_idx / blk_i;
	}

	size_t offset = 1;
	for (long j = i-1; j >= 0; j--) {
	  if (view.dim(j) == View::SplitAll) {
	    offset *= value_[j];
	  } else if (view.dim(j) > 0) {
	    offset *= view.dim(j);
	  }
	}
	idx += key_idx * offset;
      }
    }
    return idx;
#else
    // It returns row-ordered index
    size_t idx = 0;
    for (size_t i = 0; i < size_; i++) {
      if (view.dim(i) == View::SplitAll || view.dim(i) != View::SplitNone) {
	size_t key_idx = key.dim(i);
	if (view.dim(i) > 0) { // split count is specified
	  size_t blk_i = value_[i] / view.dim(i);
	  key_idx = key_idx / blk_i;
	}

	size_t offset = 1;
	for (size_t j = i+1; j < size_; j++) {
	  if (view.dim(j) == View::SplitAll) {
	    offset *= value_[j];
	  } else if (view.dim(j) > 0) {
	    offset *= view.dim(j);
	  }
	}
	idx += key_idx * offset;
      }
    }
    return idx;
#endif
  }

}

//##########################################################################//
// Method implementations common for both serial and KMR backend
//
// The following methods are copied in both two implementation headers.
//##########################################################################//
namespace kmrnext {
  using namespace std;

  template <typename DE>
  void DataStoreImpl<DE>::set(const size_t* val) {
#if VALIDATION
    if (dlist_size_ != 0) {
      throw runtime_error("DataStore is already initialized.");
    }
#endif

    dlist_size_ = 1;
    for (size_t i = 0; i < size_; i++) {
      value_[i] = val[i];
      dlist_size_ *= val[i];
    }
    DE de;
    dlist_.resize(dlist_size_, de);

    icache_.initialize(value_, dlist_size_, size_);
  }

  template <typename DE>
  void DataStoreImpl<DE>::set_dim(const size_t idx, const size_t siz) {
#if VALIDATION
    if (dlist_size_ != 0) {
      throw runtime_error("DataStore is already initialized.");
    }
#endif
    base::set_dim(idx, siz);
  }

  template <typename DE>
  void DataStoreImpl<DE>::set_from(const vector<DataStore*>& dslist) {
#if VALIDATION
    if (dslist.size() == 0) {
      throw runtime_error("There should be at least one DataStore.");
    }
    if (dlist_size_ != 0) {
      throw runtime_error("DataStore is already initialized.");
    }
    {
      // Check each DataStore in the vector
      size_t expected_dim_size = size_ - 1;
      size_t expected_dlist_size = 0;
      for (size_t i = 0; i < dslist.size(); i++) {
	DataStoreImpl<DE>* src =
	  static_cast<DataStoreImpl<DE>*>(dslist.at(i));
	if (expected_dim_size != src->size_) {
	  throw runtime_error("Dimension size of one of DataStore is wrong.");
	}
	size_t calc_dlist_size = 1;
	for (size_t j = 0; j < expected_dim_size; j++) {
	  calc_dlist_size *= src->value_[j];
	}
	if (i == 0) {
	  expected_dlist_size = calc_dlist_size;
	} else {
	  if (expected_dlist_size != calc_dlist_size) {
	    throw runtime_error("Data count of one of DataStore is wrong.");
	  }
	}
      }
    }
#endif

    size_t sizes[kMaxDimensionSize];
    sizes[0] = dslist.size();
    DataStoreImpl<DE>* ds0 = static_cast<DataStoreImpl<DE>*>(dslist.at(0));
    for (size_t i = 1; i < size_; i++) {
      sizes[i] = ds0->value_[i-1];
    }
    set(sizes);

    size_t offset = 0;
    for (size_t i = 0; i < dslist.size(); i++) {
      DataStoreImpl<DE>* src = static_cast<DataStoreImpl<DE>*>(dslist.at(i));
      for (size_t j = 0; j < src->dlist_size_; j++) {
	if (!src->dlist_[j].is_set()) {
	  continue;
	}
	dlist_[offset + j].set(src->dlist_[j].value(),
			       src->dlist_[j].size());
      }
      offset += src->dlist_size_;
    }
  }

  template <typename DE>
  void DataStoreImpl<DE>::split_to(vector<DataStore*>& dslist) {
#if VALIDATION
    if (dlist_size_ == 0) {
      throw runtime_error("Data should be set.");
    }
    if (size_ < 2) {
      throw runtime_error("DataStore can't be split.");
    }
    if (value_[0] != dslist.size()) {
      ostringstream os;
      os << "DataStore vector size should be " << value_[0]
	 << ", but " << dslist.size() << ".";
      throw runtime_error(os.str());
    }
    {
      // Check each DataStore in the vector
      size_t expected_dim_size = size_ - 1;
      for (vector<DataStore*>::iterator itr = dslist.begin();
	   itr != dslist.end(); itr++) {
	DataStoreImpl<DE>* dsi = static_cast<DataStoreImpl<DE>*>(*itr);
	if (expected_dim_size != dsi->size_) {
	  throw runtime_error("Dimension size of one of DataStore is wrong.");
	}
      }
    }
#endif

    size_t split_dims[kMaxDimensionSize];
    for (size_t i = 1; i < size_; i++) {
      split_dims[i-1] = value_[i];
    }

    size_t offset = 0;
    for (size_t i = 0; i < dslist.size(); i++) {
      DataStoreImpl<DE>* dst = static_cast<DataStoreImpl<DE>*>(dslist.at(i));
      dst->set(split_dims);
      for (size_t j = 0; j < dst->dlist_size_; j++) {
	if (!dlist_[offset + j].is_set()) {
	  continue;
	}
	dst->dlist_[j].set(dlist_[offset + j].value(),
			   dlist_[offset + j].size());
      }
      offset += dst->dlist_size_;
    }
  }

  template <typename DE>
  void DataStoreImpl<DE>::load_files(const vector<string>& files,
				     Loader<string>& loader) {
    load_array(files, loader, kmrnext_, this, value_, size_);
  }

  template <typename DE>
  void DataStoreImpl<DE>::load_integers(const vector<long>& ints,
					Loader<long>& loader) {
    load_array(ints, loader, kmrnext_, this, value_, size_);
  }

  template <typename DE>
  DataStore* DataStoreImpl<DE>::duplicate() {
    DataStoreImpl<DE>* ds = new DataStoreImpl<DE>(size_, kmrnext_);
    perform_duplicate(ds);
    return ds;
  }

  template <typename DE>
  DataElement& DataStoreImpl<DE>::data_element_at(const Key& key) {
    size_t idx = key_to_index(key);
    return dlist_[idx];
  }

}

#endif
