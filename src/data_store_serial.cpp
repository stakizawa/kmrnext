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

  void DataStore::set(const size_t *val) {
    if (dlist_size_ != 0) {
      throw runtime_error("DataStore is already initialized.");
    }

    dlist_size_ = 1;
    for (size_t i = 0; i < size_; i++) {
      value_[i] = val[i];
      dlist_size_ *= val[i];
    }
    dlist_ = new DataElement[dlist_size_];
    dlist_allocated_ = true;
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

  void DataStore::set_from(const vector<DataStore*>& dslist) {
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
	DataStore *src = dslist.at(i);
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
	    throw runtime_error("Data count one of DataStore is wrong.");
	  }
	}
      }
    }

    value_[0] = dslist.size();
    DataStore *ds0 = dslist.at(0);
    for (size_t i = 1; i < size_; i++) {
      value_[i] = ds0->value_[i-1];
    }

    dlist_size_ = 1;
    for (size_t i = 0; i < size_; i++) {
      dlist_size_ *= value_[i];
    }
    dlist_ = new DataElement[dlist_size_];
    dlist_allocated_ = true;

    size_t offset = 0;
    for (size_t i = 0; i < dslist.size(); i++) {
      DataStore *src = dslist.at(i);
      if (io_mode() == KMRNext::File) {
	src->load();
      }
      for (size_t j = 0; j < src->dlist_size_; j++) {
	dlist_[offset + j].set(src->dlist_[j].data());
      }
      if (io_mode() == KMRNext::File) {
	src->clear_cache();
      }
      offset += src->dlist_size_;
    }

    if (io_mode() == KMRNext::File) {
      store();
      clear_cache();
    }
  }

  void DataStore::split_to(vector<DataStore*>& dslist) {
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
	if (expected_dim_size != (*itr)->size_) {
	  throw runtime_error("Dimension size of one of DataStore is wrong.");
	}
      }
    }

    if (io_mode() == KMRNext::File) {
      load();
    }

    size_t split_dims[kMaxDimensionSize];
    for (size_t i = 1; i < size_; i++) {
      split_dims[i-1] = value_[i];
    }

    size_t offset = 0;
    for (size_t i = 0; i < dslist.size(); i++) {
      DataStore *dst = dslist.at(i);
      dst->set(split_dims);
      for (size_t j = 0; j < dst->dlist_size_; j++) {
	dst->dlist_[j].set(dlist_[offset + j].data());
      }
      if (io_mode() == KMRNext::File) {
	dst->store();
	dst->clear_cache();
      }
      offset += dst->dlist_size_;
    }

    if (io_mode() == KMRNext::File) {
      clear_cache();
    }
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

  DataStore* DataStore::duplicate() {
    class Copier : public Mapper {
    public:
      int operator()(DataStore *inds, DataStore *outds,
		     Key& key, vector<DataPack>& dps,
		     MapEnvironment& env)
      {
	for (std::vector<kmrnext::DataPack>:: iterator itr = dps.begin();
	     itr != dps.end(); itr++) {
	  Data d((*itr).data().value(), (*itr).data().size());
	  outds->add((*itr).key(), d);
	}
	return 0;
      }
    } copier;

    View view(size_);
    for (size_t i = 0; i < size_; i++) {
      view.set_dim(i, false);
    }
    DataStore* ds = new DataStore(size_, kmrnext_);
    ds->set(value_);
    map(copier, view, ds);
    return ds;
  }

  size_t DataStore::key_to_index(const Key& key) {
    size_t idx = 0;
    for (size_t i = 0; i < size_; i++) {
      size_t offset = 1;
      for (size_t j = i+1; j < size_; j++) {
	offset *= value_[j];
      }
      idx += key.dim(i) * offset;
    }
    return idx;
  }

  Key DataStore::index_to_key(const size_t index) {
    Key key(size_);
    size_t _index = index;
    for (size_t i = 0; i < size_; i++) {
      size_t length = 1;
      for (size_t j = i+1; j < size_; j++) {
	length *= value_[j];
      }
      key.set_dim(i, _index / length);
      _index %= length;
    }
    return key;
  }

  size_t DataStore::key_to_viewed_index(const Key& key, const View& view) {
    size_t idx = 0;
    for (size_t i = 0; i < size_; i++) {
      if (view.dim(i)) {
	size_t offset = 1;
	for (size_t j = i+1; j < size_; j++) {
	  if (view.dim(j)) {
	    offset *= value_[j];
	  }
	}
	idx += key.dim(i) * offset;
      }
    }
    return idx;
  }

  Key DataStore::key_to_viewed_key(const Key& key, const View& view) {
    size_t viewed_key_size = 0;
    for (size_t i = 0; i < size_; i++) {
      if (view.dim(i)) {
	viewed_key_size += 1;
      }
    }

    Key viewed_key(viewed_key_size);
    size_t viewed_key_idx = 0;
    for (size_t i = 0; i < size_; i++) {
      if (view.dim(i)) {
	viewed_key.set_dim(viewed_key_idx, key.dim(i));
	viewed_key_idx += 1;
      }
    }
    return viewed_key;
  }

#if 0
  // A version where dimension of the viewed key is same as the source key
  Key DataStore::key_to_viewed_key(const Key& key, const View& view) {
    Key viewed_key(size_);
    for (size_t i = 0; i < size_; i++) {
      if (view.dim(i)) {
	viewed_key.set_dim(i, key.dim(i));
      } else {
	viewed_key.set_dim(i, 0);
      }
    }
    return viewed_key;
  }
#endif

  string DataStore::filename() {
    ostringstream os;
    os << "./" << this << ".dat";
    return os.str();
  }

}
