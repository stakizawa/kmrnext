/// This file is inclured in data_store.cpp

namespace kmrnext {

  DataStore::DataStore(size_t siz, KMRNext *kn)
    : base(siz), dlist_(vector<DataElement*>()), dlist_size_(0),
      icache_(IndexCache()), map_inplace_(false), kmrnext_(kn) {}

  DataStore::~DataStore() {
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
    size_t idx = key_to_index(key);
    if (dlist_[idx] == NULL) {
      dlist_[idx] = __create_de();
    }
    if (map_inplace_) {
      dlist_[idx]->replace(data.value(), data.size());
    } else {
      dlist_[idx]->set(data.value(), data.size());
    }
  }

  DataPack DataStore::get(const Key& key) {
#if VALIDATION
    check_key_range(key);
#endif
    size_t idx = key_to_index(key);
    if (dlist_[idx] != NULL) {
      Data dat(dlist_[idx]->value(), dlist_[idx]->size());
      return DataPack(key, dat, true);
    } else {
      return DataPack(key, Data(NULL, 0), false);
    }
  }

  vector<DataPack>* DataStore::get(const View& view, const Key& key) {
#if VALIDATION
    check_view(view);
    check_key_range(key);
#endif

    size_t* blk_sizs = new size_t[size_];
    for (size_t i = 0; i < size_; i++) {
      blk_sizs[i] = (view.dim(i) > 0)? (value_[i] / view.dim(i)) : 1;
    }

    vector<DataPack> *dps = new vector<DataPack>();
    for (size_t i = 0; i < dlist_size_; i++) {
      Key tmpkey = index_to_key(i);
      bool push = true;
      for (size_t j = 0; j < size_; j++) {
	if ((view.dim(j) == View::SplitAll && key.dim(j) != tmpkey.dim(j)) ||
	    (view.dim(j) > 0 && key.dim(j) != tmpkey.dim(j) / blk_sizs[j])) {
	  push = false;
	  break;
	}
      }
      if (push && dlist_[i] != NULL) {
	Data dat(dlist_[i]->value(), dlist_[i]->size());
	dps->push_back(DataPack(tmpkey, dat, true));
      }
    }

    delete[] blk_sizs;
    return dps;
  }

  DataPack DataStore::remove(const Key& key) {
#if VALIDATION
    check_key_range(key);
#endif
    size_t idx = key_to_index(key);
    if (dlist_[idx] != NULL) {
      Data dat(dlist_[idx]->value(), dlist_[idx]->size());
      DataPack dp(key, dat, true);
      dlist_[idx]->clear();
      return dp;
    } else {
      return DataPack(key, Data(NULL, 0), false);
    }
  }

  void DataStore::map(Mapper& m, const View& view, DataStore* outds) {
#if VALIDATION
    check_map_args(view, outds);
#endif
    if (dlist_size_ == 0) {
      return;
    }

    size_t nkeys = 1;
    for (size_t i = 0; i < size_; i++) {
      if (view.dim(i) == View::SplitAll) {
	nkeys *= value_[i];
      } else if (view.dim(i) > 0) { // split count is specified
	nkeys *= view.dim(i);
      }
    }

    vector< vector<DataPack> > dpgroups(nkeys);
    for (size_t i = 0; i < dlist_size_; i++) {
      if (dlist_[i] == NULL || !dlist_[i]->is_set()) {
	continue;
      }
      Key tmpkey = index_to_key(i);
      size_t viewed_idx = key_to_viewed_index(tmpkey, view);
      vector<DataPack>& dps = dpgroups.at(viewed_idx);
      Data dat(dlist_[i]->value(), dlist_[i]->size());
      dps.push_back(DataPack(tmpkey, dat));
    }

    if (kmrnext_->profile()) {
      long dlist_count = 0;
      for (size_t i = 0; i < dpgroups.size(); i++) {
	vector<DataPack> &dps = dpgroups.at(i);
	if (dps.size() > 0) {
	  dlist_count += 1;
	}
      }
      ostringstream os;
      os << "count of data to be mapped: " << dlist_count;
      profile_out(os.str());
    }

    DataStore* _outds = outds;
    if (outds == self_ || outds == this) {
      map_inplace_ = true;
      _outds = this;
    }
    MapEnvironment env = { 0, view };
    for (size_t i = 0; i < dpgroups.size(); i++) {
      vector<DataPack>& dps = dpgroups.at(i);
      if (dps.size() > 0) {
	Key viewed_key = key_to_viewed_key(dps.at(0).key(), view);
	m(this, _outds, viewed_key, dps, env);
      }
    }
    if (outds == self_ || outds == this) {
      map_inplace_ = false;
    }
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
	result_ = os.str();
	return 0;
      }
    } dmpr(dumper);

    View view(size_);
    for (size_t i = 0; i < size_; i++) {
      view.set_dim(i, View::SplitNone);
    }

    map(dmpr, view);
    return dmpr.result_;
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
	result_ = static_cast<long>(dps.size());
	return 0;
      }
    } counter;

    View view(size_);
    for (size_t i = 0; i < size_; i++) {
      view.set_dim(i, View::SplitNone);
    }
    map(counter, view);
    return counter.result_;
  }

  DataPack SimpleFileDataStore::get(const Key& key) {
    load();
    DataPack dp = DataStore::get(key);
    clear_cache();
    return dp;
  }

  vector<DataPack>* SimpleFileDataStore::get(const View& view, const Key& key)
  {
    load();
    vector<DataPack>* dps = DataStore::get(view, key);
    clear_cache();
    return dps;
  }

  DataPack SimpleFileDataStore::remove(const Key& key) {
    load();
    DataPack dp = DataStore::remove(key);
    data_updated_ = true;
    clear_cache();
    return dp;
  }

  string SimpleFileDataStore::filename() {
    ostringstream os;
    os << "./" << this << ".dat";
    return os.str();
  }

}
