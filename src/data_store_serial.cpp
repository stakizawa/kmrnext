/// This file is inclured in data_store.cpp

namespace kmrnext {

  DataStore::DataStore(size_t siz)
    : Dimensional<size_t>(siz), dlist_(NULL), dlist_size_(0),
    dlist_allocated_(false), map_inplace_(false), kmrnext_(NULL),
    data_updated_(false), data_cached_(false) {}

  DataStore::DataStore(size_t siz, KMRNext *kn)
    : Dimensional<size_t>(siz), dlist_(NULL), dlist_size_(0),
    dlist_allocated_(false), map_inplace_(false), kmrnext_(kn),
    data_updated_(false), data_cached_(false) {}

  DataStore::~DataStore() {
    if (dlist_allocated_) {
      delete[] dlist_;
    }
    if (io_mode() == KMRNext::File) {
      string fname = filename();
      delete_file(fname);
    }
  }

  void DataStore::add(const Key& key, const Data& data) {
    check_key_range(key);
    if (!dlist_allocated_) {
      set(value_);
    }
    size_t idx = key_to_index(key);
    DataElement *de = &(dlist_[idx]);
    if (map_inplace_) {
      de->replace(&data);
    } else {
      de->set(&data);
    }
    if (io_mode() == KMRNext::File) {
      data_updated_ = true;
    }
  }

  DataPack DataStore::get(const Key& key) {
    check_key_range(key);
    if (io_mode() == KMRNext::File) {
      load();
    }
    size_t idx = key_to_index(key);
    Data* dat = dlist_[idx].data();
    DataPack dp = DataPack(key, dat, true);
    if (io_mode() == KMRNext::File) {
      clear_cache();
    }
    return dp;
  }

  vector<DataPack>* DataStore::get(const View& view, const Key& key) {
    check_view(view);
    check_key_range(key);
    if (io_mode() == KMRNext::File) {
      load();
    }

    vector<DataPack> *dps = new vector<DataPack>();
    for (size_t i = 0; i < dlist_size_; i++) {
      Key tmpkey = index_to_key(i);
      bool push = true;
      for (size_t j = 0; j < size_; j++) {
	if (view.dim(j) && key.dim(j) != tmpkey.dim(j)) {
	  push = false;
	  break;
	}
      }
      if (push) {
	dps->push_back(DataPack(tmpkey, dlist_[i].data(), true));
      }
    }

    if (io_mode() == KMRNext::File) {
      clear_cache();
    }
    return dps;
  }

  DataPack DataStore::remove(const Key& key) {
    check_key_range(key);
    if (io_mode() == KMRNext::File) {
      load();
    }

    size_t idx = key_to_index(key);
    DataPack dp(key, dlist_[idx].data(), true);
    dlist_[idx].clear();

    if (io_mode() == KMRNext::File) {
      data_updated_ = true;
      clear_cache();
    }
    return dp;
  }

  void DataStore::map(Mapper& m, const View& view, DataStore* outds) {
    check_map_args(view, outds);
    if (dlist_size_ == 0) {
      return;
    }

    if (io_mode() == KMRNext::File) {
      store();
      load();
    }

    size_t nkeys = 1;
    for (size_t i = 0; i < size_; i++) {
      if (view.dim(i)) {
	nkeys *= value_[i];
      }
    }

    vector< vector<DataPack> > dpgroups(nkeys);
    for (size_t i = 0; i < dlist_size_; i++) {
      if (dlist_[i].data() == NULL) {
	continue;
      }
      Key tmpkey = index_to_key(i);
      size_t viewed_idx = key_to_viewed_index(tmpkey, view);
      vector<DataPack>& dps = dpgroups.at(viewed_idx);
      dps.push_back(DataPack(tmpkey, dlist_[i].data()));
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

    if (io_mode() == KMRNext::File) {
      _outds->store();
      _outds->clear_cache();
      clear_cache();
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
      view.set_dim(i, false);
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
      view.set_dim(i, false);
    }
    map(counter, view);
    return counter.result_;
  }

  string DataStore::filename() {
    ostringstream os;
    os << "./" << this << ".dat";
    return os.str();
  }

}
