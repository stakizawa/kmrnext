/// This file is inclured in data_store.cpp

namespace kmrnext {

  DataStore::DataStore(size_t siz, KMRNext *kn)
    : Dimensional<size_t>(siz), dlist_(vector<DataElement*>()),
    dlist_size_(0), map_inplace_(false), kmrnext_(kn) {}

  DataStore::~DataStore() {
    if (dlist_size_ != 0) {
#ifdef _OPENMP
      #pragma omp parallel for
#endif
      for (size_t i = 0; i < dlist_size_; i++) {
	delete dlist_[i];
      }
    }
  }

  void DataStore::add(const Key& key, const Data& data) {
    check_key_range(key);
    if (dlist_size_ == 0) {
      set(value_);
    }
    size_t idx = key_to_index(key);
    DataElement *de = dlist_[idx];
    if (map_inplace_) {
      de->replace(&data);
    } else {
      de->set(&data);
    }
  }

  DataPack DataStore::get(const Key& key) {
    check_key_range(key);
    size_t idx = key_to_index(key);
    Data* dat = dlist_[idx]->data();
    return DataPack(key, dat, true);
  }

  vector<DataPack>* DataStore::get(const View& view, const Key& key) {
    check_view(view);
    check_key_range(key);

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
      if (push) {
	dps->push_back(DataPack(tmpkey, dlist_[i]->data(), true));
      }
    }

    delete[] blk_sizs;
    return dps;
  }

  DataPack DataStore::remove(const Key& key) {
    check_key_range(key);
    size_t idx = key_to_index(key);
    DataPack dp(key, dlist_[idx]->data(), true);
    dlist_[idx]->clear();
    return dp;
  }

  void DataStore::map(Mapper& m, const View& view, DataStore* outds) {
    check_map_args(view, outds);
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
      if (dlist_[i]->data() == NULL) {
	continue;
      }
      Key tmpkey = index_to_key(i);
      size_t viewed_idx = key_to_viewed_index(tmpkey, view);
      vector<DataPack>& dps = dpgroups.at(viewed_idx);
      dps.push_back(DataPack(tmpkey, dlist_[i]->data()));
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
