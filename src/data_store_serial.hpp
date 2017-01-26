#ifndef DATA_STORE_SERIAL_HPP
#define DATA_STORE_SERIAL_HPP
/// \file
/// DataStore implementation header for serial backend

#include "data_store.hpp"
#include "data_element.hpp"

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

    // The implementation is common to both Serial and KMR backend
    void load_files(const vector<string>& files, Loader<string>& loader);

    // The implementation is common to both Serial and KMR backend
    void load_integers(const vector<long>& ints, Loader<long>& loader);

    // The implementation is common to both Serial and KMR backend
    virtual DataStore* duplicate();

    string dump(DataPack::Dumper& dumper);

    long count();

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
// Method implementations specialized for the serial backend
//##########################################################################//
namespace kmrnext {
  using namespace std;

  template <typename DE>
  DataStoreImpl<DE>::DataStoreImpl(const size_t siz, KMRNext* kn)
    : base(siz), dlist_(vector<DE>()), dlist_size_(0),
      map_inplace_(false), kmrnext_(kn) {}

  template <typename DE>
  DataStoreImpl<DE>::~DataStoreImpl() {}

  template <typename DE>
  void DataStoreImpl<DE>::add(const Key& key, const Data& data) {
#if VALIDATION
    check_key_range(key);
#endif
    if (dlist_size_ == 0) {
      set(value_);
    }
    size_t idx = key_to_index(key);
    if (map_inplace_) {
      dlist_[idx].replace(data.value(), data.size());
    } else {
      dlist_[idx].set(data.value(), data.size());
    }
  }

  template <typename DE>
  DataPack DataStoreImpl<DE>::get(const Key& key) {
#if VALIDATION
    check_key_range(key);
#endif
    size_t idx = key_to_index(key);
    if (dlist_[idx].is_set()) {
      Data dat(dlist_[idx].value(), dlist_[idx].size());
      return DataPack(key, dat, true);
    } else {
      return DataPack(key, Data(NULL, 0), false);
    }
  }

  template <typename DE>
  vector<DataPack>* DataStoreImpl<DE>::get(const View& view, const Key& key) {
#if VALIDATION
    check_view(view);
    check_key_range(key);
#endif

    size_t* blk_sizs = new size_t[size_];
    for (size_t i = 0; i < size_; i++) {
      blk_sizs[i] = (view.dim(i) > 0)? (value_[i] / view.dim(i)) : 1;
    }

    vector<DataPack>* dps = new vector<DataPack>();
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
      if (push && dlist_[i].is_set()) {
	Data dat(dlist_[i].value(), dlist_[i].size());
	dps->push_back(DataPack(tmpkey, dat, true));
      }
    }

    delete[] blk_sizs;
    return dps;
  }

  template <typename DE>
  DataPack DataStoreImpl<DE>::remove(const Key& key) {
#if VALIDATION
    check_key_range(key);
#endif
    size_t idx = key_to_index(key);
    if (dlist_[idx].is_set()) {
      Data dat(dlist_[idx].value(), dlist_[idx].size());
      DataPack dp(key, dat, true);
      dlist_[idx].clear();
      return dp;
    } else {
      return DataPack(key, Data(NULL, 0), false);
    }
  }

  template <typename DE>
  void DataStoreImpl<DE>::map(Mapper& m, const View& view, DataStore* outds) {
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
      if (!dlist_[i].is_set()) {
	continue;
      }
      Key tmpkey = index_to_key(i);
      size_t viewed_idx = key_to_viewed_index(tmpkey, view);
      vector<DataPack>& dps = dpgroups.at(viewed_idx);
      Data dat(dlist_[i].value(), dlist_[i].size());
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
