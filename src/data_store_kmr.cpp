/// This file is inclured in data_store.cpp

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

namespace {
  using namespace kmrnext;

  /// Parameter for mapper_get_view
  typedef struct {
    DataStore *ds;
    vector<DataElement*>& dlist;
    vector<DataPack> *dps;
  } param_mapper_get_view;

  /// Parameter for mapper_map
  typedef struct {
    DataStore::Mapper& mapper;
    DataStore *ids;
    DataStore *ods;
    const View& view;
    DataStore::MapEnvironment& env;
    vector< vector<DataPack> >& dpgroups;
    bool map_local;
  } param_mapper_map;

  /// Parameter for mapper_collate
  typedef struct {
    DataStore *ds;
    int rank;
  } param_mapper_collate;

  /// A data set used in collate()
  class CollatePack {
  public:
    Key key_;
    DataElement *de_;
    CollatePack(Key& k, DataElement* de) : key_(k), de_(de) {}
  };

  /// It copies a key-value to vector<DataPack>.
  ///
  /// It is a KMR mapper function.  Key of the key-value should be an
  /// integer and value of the key-value should be a DataPack.  It also
  /// sets shared flags to data in the DataStore.
  int mapper_get_view(const struct kmr_kv_box kv0, const KMR_KVS *kvi,
		      KMR_KVS *kvo, void *p, const long i);

  /// It maps Datas in a DataStore.
  ///
  /// It is a KMR mapper function.  Key and value of the key-value should
  /// be an integer that represents index of Key-Data vactor.
  int mapper_map(const struct kmr_kv_box kv0, const KMR_KVS *kvi,
		 KMR_KVS *kvo, void *p, const long i);

  /// It deserializes a Data from a key-value and add the Data to
  /// the DataStore.
  ///
  /// It is a KMR mapper function.  The key of the key-value should be
  /// an integer and the value should be a byte array that represents
  /// a DataPack.
  int mapper_collate(const struct kmr_kv_box kv0, const KMR_KVS *kvi,
		     KMR_KVS *kvo, void *p, const long i);

  /// It is a wrapper of KMR map function that locally maps key-values in
  /// the input KVS.
  int map_local(KMR_KVS *kvi, KMR_KVS *kvo, void *arg, kmr_mapfn_t m);

  /// It serializes Key, Data and their attributes to a byte array
  /// for collate().
  void serialize_collate(const CollatePack* dep, char** buf, size_t* buf_siz);

  /// It deserializes a Key, Data and their attribute from a byte array
  /// for collate().
  void deserialize_collate(const char* buf, size_t buf_siz, DataStore* ds,
			   int owner);

  /// It calculates send target rank.
  inline size_t calc_send_target(size_t index, size_t count, size_t nprocs);
}

namespace kmrnext {

  DataStore::DataStore(size_t siz, KMRNext *kn)
    : base(siz), dlist_(vector<DataElement*>()),
      dlist_size_(0), icache_(IndexCache()), map_inplace_(false), kmrnext_(kn),
      parallel_(false), split_(NULL), collated_(false) {}

  DataStore::~DataStore() {
    if (split_ != NULL) {
      delete split_;
    }
    if (dlist_size_ != 0) {
#ifdef _OPENMP
      #pragma omp parallel for schedule(static, OMP_FOR_CHUNK_SIZE)
#endif
      for (size_t i = 0; i < dlist_size_; i++) {
	if (dlist_[i] != NULL) {
	  delete dlist_[i];
	}
      }
    }
  }

  void DataStore::add(const Key& key, const Data& data) {
#if VALIDATION
    check_key_range(key);
#endif
    if (dlist_size_ == 0) {
      set(value_);
    }
    if (parallel_ || kmrnext_->rank() == 0) {
      size_t idx = key_to_index(key);
      if (dlist_[idx] == NULL) {
	dlist_[idx] = __create_de();
      }
      try {
	if (map_inplace_) {
	  dlist_[idx]->replace(&data);
	} else {
	  dlist_[idx]->set(&data);
	}
      }
      catch (runtime_error& e) {
	cerr << "Failed to add a data to DataStore" << to_string()
	     << " at Key" << key.to_string()
	     << "on Rank" << kmrnext_->rank() << "." << endl;
	throw e;
      }
      dlist_[idx]->set_owner(kmrnext_->rank());
    }
  }

  DataPack DataStore::get(const Key& key) {
#if VALIDATION
    check_key_range(key);
#endif

    size_t idx = key_to_index(key);
    if (dlist_[idx] == NULL) {
      dlist_[idx] = __create_de();
    }
    if (dlist_[idx]->is_shared()) {
      return DataPack(key, dlist_[idx]->data(), true);
    }

    KMR_KVS *snd = kmr_create_kvs(kmrnext_->kmr(),
				  KMR_KV_INTEGER, KMR_KV_OPAQUE);
    KMR_KVS *rcv = kmr_create_kvs(kmrnext_->kmr(),
				  KMR_KV_INTEGER, KMR_KV_OPAQUE);
    if (dlist_[idx]->is_set()) {
      Data *d = dlist_[idx]->data();
      size_t snd_siz = sizeof(int) + d->size();
      char *snd_buf = static_cast<char*>(malloc(sizeof(char) * snd_siz));
      int owner = kmrnext_->rank();
      memcpy(snd_buf, &owner, sizeof(int));
      memcpy(snd_buf + sizeof(int), d->value(), d->size());
      struct kmr_kv_box kv;
      kv.klen = (int)sizeof(long);
      kv.vlen = (int)snd_siz;
      kv.k.i  = idx;
      kv.v.p  = snd_buf;
      kmr_add_kv(snd, kv);
      free(snd_buf);
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
      if (!dlist_[idx]->is_set()) {
	// Copy data for future access.
	int owner;
	memcpy(&owner, kv.v.p, sizeof(int));
	Data rcvdat(static_cast<void*>((const_cast<char*>(kv.v.p) +
					sizeof(int))),
		    kv.vlen - sizeof(int));
	dlist_[idx]->set(&rcvdat);
	dlist_[idx]->set_owner(owner);
      }
      dlist_[idx]->shared();
    }
    kmr_free_kvs(rcv);

    return DataPack(key, dlist_[idx]->data(), true);
  }

  vector<DataPack>* DataStore::get(const View& view, const Key& key) {
#if VALIDATION
    check_view(view);
    check_key_range(key);
#endif

    // Allocate all DataElement beforehand
#ifdef _OPENMP
    #pragma omp parallel for schedule(static, OMP_FOR_CHUNK_SIZE)
#endif
    for (size_t i = 0; i < dlist_size_; i++) {
      if (dlist_[i] == NULL) {
	dlist_[i] = __create_de();
      }
    }

    vector<DataPack> *dps = new vector<DataPack>();

    size_t* blk_sizs = new size_t[size_];
    for (size_t i = 0; i < size_; i++) {
      blk_sizs[i] = (view.dim(i) > 0)? (value_[i] / view.dim(i)) : 1;
    }

    KMR_KVS *snd = kmr_create_kvs(kmrnext_->kmr(),
				  KMR_KV_INTEGER, KMR_KV_OPAQUE);
    KMR_KVS *rcv = kmr_create_kvs(kmrnext_->kmr(),
				  KMR_KV_INTEGER, KMR_KV_OPAQUE);
#ifdef _OPENMP
    #pragma omp parallel for schedule(static, OMP_FOR_CHUNK_SIZE)
#endif
    for (size_t i = 0; i < dlist_size_; i++) {
      if (!dlist_[i]->is_set()) {
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
	if (dlist_[i]->is_shared()) {
	  // the data is already replicated
#ifdef _OPENMP
          #pragma omp critical
#endif
	  dps->push_back(DataPack(tmpkey, dlist_[i]->data(), true));
	  continue;
	} else {
	  // replicate the data
	  Data *d = dlist_[i]->data();
	  size_t snd_siz = sizeof(int) + d->size();
	  char *snd_buf = static_cast<char*>(malloc(sizeof(char) * snd_siz));
	  int owner = kmrnext_->rank();
	  memcpy(snd_buf, &owner, sizeof(int));
	  memcpy(snd_buf + sizeof(int), d->value(), d->size());
	  struct kmr_kv_box kv;
	  kv.klen = static_cast<int>(sizeof(size_t));
	  kv.vlen = static_cast<int>(snd_siz);
	  kv.k.i  = i;
	  kv.v.p  = snd_buf;
	  kmr_add_kv(snd, kv);
	  free(snd_buf);
	}
      }
    }
    kmr_add_kv_done(snd);
    kmr_replicate(snd, rcv, kmr_noopt);
    param_mapper_get_view param = { this, dlist_, dps };
    kmr_map(rcv, NULL, (void*)&param, kmr_noopt, mapper_get_view);

    delete[] blk_sizs;
    return dps;
  }

  DataPack DataStore::remove(const Key& key) {
    DataPack dp0 = get(key);
    if (dp0.data().value() != NULL) {
      size_t idx = key_to_index(key);
      dlist_[idx]->clear();
    }
    return dp0;
  }

  void DataStore::map(Mapper& m, const View& view, DataStore* outds) {
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
    #pragma omp parallel for schedule(static, OMP_FOR_CHUNK_SIZE)
#endif
    for (size_t i = 0; i < dlist_size_; i++) {
      if (dlist_[i] == NULL ||
	  dlist_[i]->data() == NULL ||
	  (dlist_[i]->is_shared() && dlist_[i]->owner() != kmrnext_->rank())) {
	continue;
      }
      Key tmpkey = index_to_key(i);
      size_t viewed_idx = key_to_viewed_index(tmpkey, view);
      vector<DataPack>& dps = dpgroups.at(viewed_idx);
#ifdef _OPENMP
      #pragma omp critical
#endif
      dps.push_back(DataPack(tmpkey, dlist_[i]->data()));
    }

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

    KMR_KVS *ikvs = kmr_create_kvs(kmrnext_->kmr(),
				   KMR_KV_INTEGER, KMR_KV_INTEGER);
#ifdef _OPENMP
    #pragma omp parallel for
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

    DataStore *_outds = outds;
    if (outds == self_ || outds == this) {
      map_inplace_ = true;
      _outds = this;
    }
    _outds->parallel_ = true;
    MapEnvironment env = { kmrnext_->rank(), view, split,
			   MPI_COMM_NULL };
    param_mapper_map param = { m, this, _outds, view, env, dpgroups,
			       is_local };
    if (is_local) {
      double timestamp[2] = {0, 0};
      long task_count;
      kmr_local_element_count(ikvs, &task_count);
      if (kmrnext_->profile()) {
	timestamp[0] = gettime(kmrnext_->kmr()->comm);
      }

      map_local(ikvs, NULL, (void*)&param, mapper_map);

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
				  kmr_noopt, kmrnext_->rank(), mapper_map);
    }
    _outds->parallel_ = false;
    if (outds == self_ || outds == this) {
      map_inplace_ = false;
      // unshare all data
#ifdef _OPENMP
      #pragma omp parallel for schedule(static, OMP_FOR_CHUNK_SIZE)
#endif
      for (size_t i = 0; i < dlist_size_; i++) {
	if (dlist_[i] == NULL) {
	  continue;
	}
	if (dlist_[i]->is_shared() && dlist_[i]->owner() != kmrnext_->rank()) {
	  dlist_[i]->clear();
	  continue;
	}
	dlist_[i]->unshared();
      }
    }
  }

  void DataStore::set_split(const View& view) {
#if VALIDATION
    check_view(view);
#endif
    if (split_ == NULL) {
      split_ = new View(view.size());
    }
    for (size_t i = 0; i < view.size(); i++) {
      split_->set_dim(i, view.dim(i));
    }
  }

  View DataStore::get_split() {
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

  string DataStore::dump(DataPack::Dumper& dumper) {
    class WrappedDumper : public Mapper {
    public:
      string result_;
      DataPack::Dumper& dumper_;

      WrappedDumper(DataPack::Dumper& dmpr) : dumper_(dmpr) {}
      int operator()(DataStore *inds, DataStore *outds,
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
	char *local_cstr = const_cast<char*>(dumped.c_str());
	int nprocs;
	MPI_Comm_size(env.mpi_comm, &nprocs);
	int *local_lens = static_cast<int*>(malloc(sizeof(int) * nprocs));
	MPI_Allgather(&local_len, 1, MPI_INT, local_lens, 1, MPI_INT,
		      env.mpi_comm);
	int total_len = 0;
	int *displs = static_cast<int*>(malloc(sizeof(int) * nprocs));
	for (int i = 0; i < nprocs; i++) {
	  displs[i] = total_len;
	  total_len += local_lens[i];
	}
	total_len += 1; // +1 for '\0'
	char *total_cstr = static_cast<char*>(malloc(sizeof(char) * total_len));
	MPI_Allgatherv(local_cstr, local_len, MPI_CHAR,
		       total_cstr, local_lens, displs, MPI_CHAR, env.mpi_comm);
	total_cstr[total_len - 1] = '\0';

	result_.append(total_cstr);
	free(total_cstr);
	free(displs);
	free(local_lens);
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
    char *result_cstr = static_cast<char*>(malloc(sizeof(char) * length));
    if (kmrnext_->rank() == master) {
      memcpy(result_cstr, dmpr.result_.c_str(), sizeof(char) * length);
    }
    MPI_Bcast(result_cstr, length, MPI_CHAR, master, kmrnext_->kmr()->comm);
    string result(result_cstr);
    free(result_cstr);
    return result;
  }

  long DataStore::count() {
    class Counter : public Mapper {
    public:
      long result_;
      Counter() : result_(0) {}
      int operator()(DataStore *inds, DataStore *outds,
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

  size_t DataStore::key_to_split_index(const Key& key, const View& view) {
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

  void DataStore::collate() {
    double time_collate_start = 0, time_collate_finish = 0;
    if (kmrnext_->profile()) {
      time_collate_start = gettime(kmrnext_->kmr()->comm);
    }

    // Split pattern set to this DataStore
    View split = get_split();
    // Count of viewed data
    size_t ndata;
    // Indices of viewed Keys whose related data should be reside
    // in this process
    size_t *indices;
    // size of indices
    size_t indices_cnt;

    // Check if collate is necessary.
    {
      ndata = 1;
      // Element count of each data
      size_t each_cnt = 1;
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
	indices_cnt =
	  (kmrnext_->rank() < static_cast<int>(rem_cnt))?
	  avg_cnt + 1 : avg_cnt;
	if (indices_cnt == 0) {
	  indices = NULL;
	} else {
	  indices = new size_t[indices_cnt];
	  size_t start = avg_cnt * static_cast<size_t>(kmrnext_->rank())
	    + ((kmrnext_->rank() < static_cast<int>(rem_cnt))?
	       kmrnext_->rank() : rem_cnt);
	  for (size_t i = 0; i < indices_cnt; i++) {
	    indices[i] = start + i;
	  }
	}
      }

      // It holds the DataPacks this process should hold
      vector< vector<DataPack> > dpgroups(indices_cnt);
      if (indices_cnt != 0) {
	size_t range_start = indices[0];
	size_t range_end   = indices[indices_cnt - 1];
#ifdef _OPENMP
        #pragma omp parallel for schedule(static, OMP_FOR_CHUNK_SIZE)
#endif
	for (size_t i = 0; i < dlist_size_; i++) {
	  if (dlist_[i] == NULL ||
	      dlist_[i]->data() == NULL ||
	      (dlist_[i]->is_shared() &&
	       dlist_[i]->owner() != kmrnext_->rank())) {
	    continue;
	  }
	  Key tmpkey = index_to_key(i);
	  size_t viewed_idx = key_to_split_index(tmpkey, split);
	  if (range_start <= viewed_idx && viewed_idx <= range_end) {
	    size_t idx = viewed_idx - range_start;
	    vector<DataPack>& dps = dpgroups.at(idx);
#ifdef _OPENMP
            #pragma omp critical
#endif
	    dps.push_back(DataPack(tmpkey, dlist_[i]->data()));
	  }
	}
      }

      int local_need_collate = 0;
#ifdef _OPENMP
      #pragma omp parallel for
#endif
      for (size_t i = 0; i < dpgroups.size(); i++) {
	vector<DataPack>& dps = dpgroups.at(i);
	if (dps.size() != each_cnt) {
	  local_need_collate = 1;
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
	if (indices != NULL) {
	  delete[] indices;
	}
	return;
      } else {
	collated_ = true;
      }
    }

    // Collate from here

    // Search data should be sent
    vector< vector<CollatePack> > cpgroups(ndata);
    {
      size_t range_start = (indices_cnt > 0)? indices[0] : 0;
      size_t range_end   = (indices_cnt > 0)? indices[indices_cnt - 1] : 0;
#ifdef _OPENMP
      #pragma omp parallel for schedule(static, OMP_FOR_CHUNK_SIZE)
#endif
      for (size_t i = 0; i < dlist_size_; i++) {
	if (dlist_[i] == NULL || !dlist_[i]->is_set()) {
	  continue;
	}
	if (dlist_[i]->is_shared() && dlist_[i]->owner() != kmrnext_->rank()) {
	  dlist_[i]->clear();
	  continue;
	}
	dlist_[i]->unshared();
	Key tmpkey = index_to_key(i);
	size_t viewed_idx = key_to_split_index(tmpkey, split);
	if (indices_cnt > 0 &&
	    !(range_start <= viewed_idx && viewed_idx <= range_end)) {
	  vector<CollatePack>& cps = cpgroups.at(viewed_idx);
#ifdef _OPENMP
          #pragma omp critical
#endif
	  cps.push_back(CollatePack(tmpkey, dlist_[i]));
	}
      }
    }

    KMR_KVS *kvs0 = kmr_create_kvs(kmrnext_->kmr(),
				   KMR_KV_INTEGER, KMR_KV_OPAQUE);
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (size_t i = 0; i < ndata; i++) {
      vector<CollatePack> &cps = cpgroups.at(i);
      if (cps.size() > 0) {
	for (vector<CollatePack>::iterator itr = cps.begin(); itr != cps.end();
	     itr++) {
	  char *buf;
	  size_t buf_siz;
	  serialize_collate(&(*itr), &buf, &buf_siz);
	  struct kmr_kv_box kv;
	  kv.klen = static_cast<int>(sizeof(size_t));
	  kv.vlen = static_cast<int>(buf_siz);
	  kv.k.i  = calc_send_target(i, ndata,
				     static_cast<size_t>(kmrnext_->nprocs()));
	  kv.v.p  = buf;
	  kmr_add_kv(kvs0, kv);
	  free(buf);
	  (*itr).de_->clear();
	}
      }
    }
    kmr_add_kv_done(kvs0);

    KMR_KVS *kvs1 = kmr_create_kvs(kmrnext_->kmr(),
				   KMR_KV_INTEGER, KMR_KV_OPAQUE);
    struct kmr_option kmr_shflopt = kmr_noopt;
    kmr_shflopt.key_as_rank = 1;
    kmr_shuffle(kvs0, kvs1, kmr_shflopt);
    parallel_ = true;
    map_inplace_ = true;
    param_mapper_collate p0 = { this, kmrnext_->rank() };
    kmr_map(kvs1, NULL, &p0, kmr_noopt, mapper_collate);
    map_inplace_ = false;
    parallel_ = false;

    if (indices != NULL) {
      delete[] indices;
    }

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

  bool DataStore::collated() {
    return collated_;
  }

  void DataStore::load_local_data(Loader<long>& loader) {
    // Create a 1D DataStore that stores ranks
    DataStore* ds0 = new DataStore(1, kmrnext_);
    ds0->set_dim(0, kmrnext_->nprocs());
    Key key(1);
    for (long i = 0; i < kmrnext_->nprocs(); i++) {
      key.set_dim(0, i);
      Data dat(&i, sizeof(long));
      ds0->add(key, dat);
    }

    // Use Split <T> so that each process loads an array element.
    View split_ds0(1);
    split_ds0.set_dim(0, View::SplitAll);
    ds0->set_split(split_ds0);

    // Define a mapper for the loader
    class WrappedLoader : public DataStore::Mapper {
    public:
      DataStore::Loader<long>& loader_;
      KMRNext* kmrnext_;

      WrappedLoader(DataStore::Loader<long>& ldr, KMRNext* next)
	: loader_(ldr), kmrnext_(next) {}
      int operator()(DataStore *inds, DataStore *outds,
		     Key& k, vector<DataPack>& dps,
		     DataStore::MapEnvironment& env)
      {
#ifdef DEBUG
	if (dps.size() != 1) {
	  throw runtime_error("System error: Unexpected shuffle behavior.");
	}
#endif
	long *dsval = static_cast<long*>(dps[0].data().value());
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
    ds0->map(wloader, v, this);
    delete ds0;

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

  DataPack SimpleFileDataStore::get(const Key& key) {
#if VALIDATION
    check_key_range(key);
#endif
    load();
    size_t idx = key_to_index(key);
    if (dlist_[idx] == NULL) {
      dlist_[idx] = __create_de();
    }
    bool unset = !(dlist_[idx]->is_set());
    DataPack dp = DataStore::get(key);
    if (unset) {
      // It assume that unset DataElement is set by the above get()
      data_updated_ = true;
      store();
    }
    clear_cache();
    return dp;
  }

  vector<DataPack>* SimpleFileDataStore::get(const View& view, const Key& key)
  {
#if VALIDATION
    check_key_range(key);
#endif
    load();
    vector<DataPack>* dps = DataStore::get(view, key);
    data_updated_ = true; // Though sometimes not updated.
    store();
    clear_cache();
    return dps;
  }

  DataPack SimpleFileDataStore::remove(const Key& key) {
    load();
    DataPack dp0 = DataStore::get(key);
    if (dp0.data().value() != NULL) {
      size_t idx = key_to_index(key);
      dlist_[idx]->clear();
      data_updated_ = true;
    }
    clear_cache();
    return dp0;
  }

  void SimpleFileDataStore::collate() {
    // It does not clear cache as callete() is usually called followed by map()
    load();
    DataStore::collate();
    // As data is exchanged between processes, the saved file is stale.
    string fname = filename();
    delete_file(fname);
    store();
  }

  string SimpleFileDataStore::filename() {
    ostringstream os;
    os << "./" << this << "_" << kmrnext_->rank() << ".dat";
    return os.str();
  }

}

namespace {
  int mapper_get_view(const struct kmr_kv_box kv0, const KMR_KVS *kvi,
		      KMR_KVS *kvo, void *p, const long i) {
    param_mapper_get_view *param = (param_mapper_get_view *)p;
    size_t idx = kv0.k.i;
    Key key = param->ds->index_to_key(idx);
    if (!param->dlist[idx]->is_set()) {
      int owner;
      memcpy(&owner, kv0.v.p, sizeof(int));
      Data data(static_cast<void*>(const_cast<char*>(kv0.v.p) + sizeof(int)),
		kv0.vlen - sizeof(int));
      param->dlist[idx]->set(&data);
      param->dlist[idx]->set_owner(owner);
    }
    param->dlist[idx]->shared();
#ifdef _OPENMP
    // As a KMR map function is run in parallel by OpenMP, shared resources
    // should be modified in critical regions.
    #pragma omp critical
#endif
    param->dps->push_back(DataPack(key, param->dlist[idx]->data(), true));
    return MPI_SUCCESS;
  }

  int mapper_map(const struct kmr_kv_box kv0, const KMR_KVS *kvi,
		 KMR_KVS *kvo, void *p, const long i) {
    param_mapper_map *param = static_cast<param_mapper_map*>(p);
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

  int map_local(KMR_KVS *kvi, KMR_KVS *kvo, void *arg, kmr_mapfn_t m) {
    param_mapper_map *param = static_cast<param_mapper_map*>(arg);
    MPI_Comm self_comm;
    MPI_Comm_split(kvi->c.mr->comm, kvi->c.mr->rank, kvi->c.mr->rank,
		   &self_comm);
    param->env.mpi_comm = self_comm;
    struct kmr_option kmr_nothrd = kmr_noopt;
    kmr_nothrd.nothreading = 1;
    kmr_map(kvi, kvo, arg, kmr_nothrd, mapper_map);
    MPI_Comm_free(&self_comm);
    return MPI_SUCCESS;
  }

  int mapper_collate(const struct kmr_kv_box kv0, const KMR_KVS *kvi,
		     KMR_KVS *kvo, void *p, const long i) {
    param_mapper_collate *param = static_cast<param_mapper_collate*>(p);
    deserialize_collate(kv0.v.p, kv0.vlen, param->ds, param->rank);
    return MPI_SUCCESS;
  }

  void serialize_collate(const CollatePack* cp, char** buf, size_t* buf_siz) {
    // Memory layout
    // 1. [key]  size    : sizeof(size_t)
    // 2. [key]  content : sizeof(size_t) * count_of_key
    // 3. [data] size    : sizeof(size_t)
    // 4. [data] content : size of 3
    // 5. [data] owner   : sizeof(int)
    //
    // Two parameters of Data, value_allocated_ and shared_ is not serialized.

    // calculate size
    Key k = cp->key_;
    Data d = cp->de_->data();
    *buf_siz  = sizeof(size_t);            // 1. [key]  size
    *buf_siz += sizeof(size_t) * k.size(); // 2. [key]  content
    *buf_siz += sizeof(size_t);            // 3. [data] size
    *buf_siz += d.size();                  // 4. [data] content
    *buf_siz += sizeof(int);               // 5. [data] owner

    // set data to buf
    *buf = static_cast<char*>(calloc(*buf_siz, sizeof(char)));
    // key
    char *p = *buf;
    *reinterpret_cast<size_t*>(p) = k.size();
    p += sizeof(size_t);
    for (size_t i = 0; i < k.size(); i++) {
      *reinterpret_cast<size_t*>(p) = k.dim(i);
      p += sizeof(size_t);
    }
    // data
    *reinterpret_cast<size_t*>(p) = d.size();
    p += sizeof(size_t);
    memcpy(p, d.value(), d.size());
    p += d.size();
    *reinterpret_cast<int*>(p) = cp->de_->owner();
  }

  void deserialize_collate(const char* buf, size_t buf_siz, DataStore* ds,
			   int owner) {
    // Memory layout
    // 1. [key]  size    : sizeof(size_t)
    // 2. [key]  content : sizeof(size_t) * count_of_key
    // 3. [data] size    : sizeof(size_t)
    // 4. [data] content : size of 3
    // 5. [data] owner   : sizeof(int)

    // key
    char *p = const_cast<char*>(buf);
    size_t key_siz = *reinterpret_cast<size_t*>(p);
    p += sizeof(size_t);
    size_t *key_ary = new size_t[key_siz];
    for (size_t i = 0; i < key_siz; i++) {
      key_ary[i] = *reinterpret_cast<size_t*>(p);
      p += sizeof(size_t);
    }
    Key k(key_siz);
    k.set(key_ary);
    delete[] key_ary;
    // data
    size_t dat_siz = *reinterpret_cast<size_t*>(p);
    p += sizeof(size_t);
    Data d(static_cast<void*>(p), dat_siz);
#if 0
    // owner in serialized data is omitted
    p += dat_siz;
    int owner = *reinterpret_cast<int*>(p);
#endif

    DataElement* de = ds->data_element_at(k);
    de->set(&d);
    de->set_owner(owner);
  }

  inline size_t calc_send_target(size_t index, size_t count, size_t nprocs) {
    size_t avg = count / nprocs;
    size_t rem = count % nprocs;
    size_t idx_max = 0;
    size_t result;
#ifdef DEBUG
    size_t fset = false;
#endif
    for (size_t i = 0; i < nprocs; i++) {
      idx_max += (i < rem)? avg + 1 : avg;
      if (index < idx_max) {
	result = i;
#ifdef DEBUG
	fset = true;
#endif
	break;
      }
    }
#ifdef DEBUG
    if (!fset) {
      throw runtime_error("Logical error at calc_send_target().");
    }
#endif
    return result;
  }

}
