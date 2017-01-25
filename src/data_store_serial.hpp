#ifndef DATA_STORE_SERIAL_HPP
#define DATA_STORE_SERIAL_HPP
/// \file
/// Store implementation header for serial backend

#include <fstream>
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

    virtual ~DataStoreImpl() {}

    void set(const size_t* val) {
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

    void set_dim(const size_t idx, const size_t siz) {
#if VALIDATION
      if (dlist_size_ != 0) {
	throw runtime_error("DataStore is already initialized.");
      }
#endif
      base::set_dim(idx, siz);
    }

    virtual void add(const Key& key, const Data& data) {
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

    virtual DataPack get(const Key& key) {
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

    virtual vector<DataPack>* get(const View& view, const Key& key) {
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

    virtual DataPack remove(const Key& key) {
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

    // The implementation is common to both Serial and KMR backend
    virtual void set_from(const vector<DataStore*>& dslist) {
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

    // The implementation is common to both Serial and KMR backend
    virtual void split_to(vector<DataStore*>& dslist) {
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

    virtual void map(Mapper& m, const View& view, DataStore* outds=self_) {
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

// TODO delete
#ifdef BACKEND_KMR
    void set_split(const View& split) {}

    View get_split() {}

    virtual void collate() {}

    bool collated() { return false; }
#endif

    // The implementation is common to both Serial and KMR backend
    void load_files(const vector<string>& files, Loader<string>& loader) {
      load_array(files, loader, kmrnext_, this, value_, size_);
    }

    // The implementation is common to both Serial and KMR backend
    void load_integers(const vector<long>& ints, Loader<long>& loader) {
      load_array(ints, loader, kmrnext_, this, value_, size_);
    }

// TODO delete
#ifdef BACKEND_KMR
    void load_parallel(Loader<long>& loader) {}
#endif

    virtual DataStore* duplicate() {
      DataStoreImpl<DE>* ds = new DataStoreImpl<DE>(size_, kmrnext_);
      perform_duplicate(ds);
      return ds;
    }

    string dump(DataPack::Dumper& dumper) {
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

    long count() {
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

#ifdef BACKEND_KMR
    size_t key_to_split_index(const Key& key, const View& view) { return 0; }
#endif

    DataElement& data_element_at(const Key& key) {
      size_t idx = key_to_index(key);
      return dlist_[idx];
    }

  protected:
    // Stored DataElements
    vector<DE> dlist_;
    // Size of dlist_
    size_t dlist_size_;
    // True if the input and output DataStore of map function is same
    bool map_inplace_;
    // A KMRNext object that stores execution status
    KMRNext* kmrnext_;
#ifdef BACKEND_KMR
    // True if the DataStore should be processed in parallel
    bool parallel_;
    // Split pattern of the DataStore, that defines data distribution
    View* split_;
    // Set to be true if the last call of map() or collate() actually
    // performed collate.
    bool collated_;
#endif

    // It creates a new DataStore.
    //
    // Size of each dimension should be set by DataStore.set() function
    // before using this DataStore.
    //
    // \param[in] size number of dimensions of this DataStore.
    // \param[in] kn   an instance of KMRNext context
    DataStoreImpl(const size_t siz, KMRNext* kn)
      : base(siz), dlist_(vector<DE>()), dlist_size_(0),
	map_inplace_(false), kmrnext_(kn) {}

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
    virtual ~SimpleFileDataStore() {
      string fname = filename();
      delete_file(fname);
    }

    // The implementation is common to both Serial and KMR backend
    void add(const Key& key, const Data& data) {
      base::add(key, data);
      data_updated_ = true;
    }

    DataPack get(const Key& key) {
      load();
      DataPack dp = base::get(key);
      clear_cache();
      return dp;
    }

    vector<DataPack>* get(const View& view, const Key& key) {
      load();
      vector<DataPack>* dps = base::get(view, key);
      clear_cache();
      return dps;
    }

    DataPack remove(const Key& key) {
      load();
      DataPack dp = base::remove(key);
      data_updated_ = true;
      clear_cache();
      return dp;
    }

    // The implementation is common to both Serial and KMR backend
    void set_from(const vector<DataStore*>& dslist) {
#ifdef _OPENMP
      #pragma omp parallel for
#endif
      for (size_t i = 0; i < dslist.size(); i++) {
	dynamic_cast<SimpleFileDataStore*>(dslist[i])->load();
      }
      base::set_from(dslist);
      store();
      clear_cache();
#ifdef _OPENMP
      #pragma omp parallel for
#endif
      for (size_t i = 0; i < dslist.size(); i++) {
	dynamic_cast<SimpleFileDataStore*>(dslist[i])->clear_cache();
      }
    }

    // The implementation is common to both Serial and KMR backend
    void split_to(vector<DataStore*>& dslist) {
      load();
      base::split_to(dslist);
      for (size_t i = 0; i < dslist.size(); i++) {
	SimpleFileDataStore* ds =
	  dynamic_cast<SimpleFileDataStore*>(dslist[i]);
	ds->store();
	ds->clear_cache();
      }
      clear_cache();
    }

    // The implementation is common to both Serial and KMR backend
    void map(Mapper& m, const View& view, DataStore* outds=self_) {
      store();
      load();
      base::map(m, view, outds);
      SimpleFileDataStore* _outds = dynamic_cast<SimpleFileDataStore*>(outds);
      if (outds == self_ || outds == this) {
	_outds = this;
	data_updated_ = true;
      }
      _outds->store();
      _outds->clear_cache();
      clear_cache();
    }

    // The implementation is common to both Serial and KMR backend
    DataStore* duplicate()  {
      SimpleFileDataStore* ds = new SimpleFileDataStore(size_, kmrnext_);
      perform_duplicate(ds);
      return ds;
    }

    // TODO delete
#ifdef BACKEND_KMR
    virtual void collate();
#endif

  private:
    // True if the data in DataStore is updated, but not written to a file
    bool data_updated_;
    // True if the Data in a file is cached in memory
    bool data_cached_;

    // It returns name of file that stores data elements of the DataStore.
    string filename() {
      ostringstream os;
      os << "./" << this << ".dat";
      return os.str();
    }

    // It writes data elements in the DataStore to a file.
    // The implementation is common to both Serial and KMR backend
    //
    // \return   true if Data are stored in a file
    bool store() {
      string fname = filename();
      if (file_exist(fname)) {
	if (data_updated_) {
	  load();
	  delete_file(fname);
	} else {
	  return false;
	}
      }

      ofstream fout;
      fout.open(fname.c_str(), ios::out|ios::binary);
      size_t write_offset = 0;

      for (size_t i = 0; i < dlist_size_; i++) {
	if (!dlist_[i].is_set()) {
	  continue;
	}
	size_t d_siz = dlist_[i].size();
	if (d_siz == 0) {
	  continue;
	}
	fout.write(static_cast<char*>(dlist_[i].value()),
		   static_cast<streamsize>(d_siz));
	dlist_[i].written(write_offset, d_siz);
	write_offset += d_siz;
      }
      fout.flush();
      fout.close();
      data_updated_ = false;
      data_cached_ = true;
      return true;
    }

    // It reads data elements in a file and then loads to the DataStore.
    // The implementation is common to both Serial and KMR backend
    //
    // \return   true if Data are loaded from a file
    bool load() {
      string fname = filename();
      if (data_cached_ && !file_exist(fname)) {
	throw runtime_error("File is not found.");
      }
      if (data_cached_ || (!data_cached_ && !file_exist(fname))) {
	return false;
      }

      size_t file_siz = file_size(fname);
      if (file_siz == 0) {
	return false;
      }
      char* buf = new char[file_siz];
      ifstream fin;
      fin.open(fname.c_str(), ios::in|ios::binary);
      fin.read(buf, static_cast<streamsize>(file_siz));
      fin.close();

#ifdef _OPENMP
      #pragma omp parallel for schedule(static, OMP_FOR_CHUNK_SIZE)
#endif
      for (size_t i = 0; i < dlist_size_; i++) {
	if (!dlist_[i].is_set()) {
	  continue;
	}
	dlist_[i].restore(buf);
      }
      delete[] buf;
      data_cached_ = true;
      return true;
    }

    // It clears memory cache.
    // The implementation is common to both Serial and KMR backend
    void clear_cache() {
      if (!data_cached_) {
	return;
      }
#ifdef _OPENMP
      #pragma omp parallel for schedule(static, OMP_FOR_CHUNK_SIZE)
#endif
      for (size_t i = 0; i < dlist_size_; i++) {
	if (!dlist_[i].is_set()) {
	  continue;
	}
	dlist_[i].clear_cache();
      }
      data_cached_ = false;
    }

  };

}

#endif
