/// \file
/// Implementation of DataStore class
///
/// Classes and methods defined in this file is common for both serial and
/// KMR backend.

#include "../config.hpp"
#include "kmrnext.hpp"
#include "data_store.hpp"
#include <fstream>

#ifdef BACKEND_SERIAL
#include "data_store_serial.hpp"
#elif defined BACKEND_KMR
#include "data_store_kmr.hpp"
#endif

namespace {
  using namespace kmrnext;

  ////////////////////////////////////////////////////////////////////////////
  // A data loader class that sets long integer 0 to all data elements.
  ////////////////////////////////////////////////////////////////////////////
  class Zeroizer : public DataStore::Loader<long> {
  public:
    int operator()(DataStore* ds, const long& num)
    {
      if (ds->size() == 1) {
	Key key(1);
	key.set_dim(0, num);
	long val = 0;
	Data data(&val, sizeof(long));
	ds->add(key, data);
      } else {
#ifdef _OPENMP
        #pragma omp parallel
#endif
	{
	  Key key(ds->size());
	  key.set_dim(0, num);
#ifdef _OPENMP
          #pragma omp for
#endif
	  for (size_t i = 0; i < ds->dim(1); i++) {
	    key.set_dim(1, i);
	    set_data(ds, key, 2, ds->size());
	  }
	}
      }
      return 0;
    }
  private:
    void set_data(DataStore* ds, Key& key, size_t cur_depth, size_t max_depth)
    {
      if (cur_depth < max_depth) {
	for (size_t i = 0; i < ds->dim(cur_depth); i++) {
	  key.set_dim(cur_depth, i);
	  set_data(ds, key, cur_depth + 1, max_depth);
	}
      } else {
	long val = 0;
	Data data(&val, sizeof(long));
	ds->add(key, data);
      }
    }
  };

}

namespace kmrnext {
  ////////////////////////////////////////////////////////////////////////////
  // Method implementation for class DataStore
  ////////////////////////////////////////////////////////////////////////////

  DataStore* DataStore::self_;

  void DataStore::initialize(KMRNext* next) {
    if (DataStore::self_ == NULL) {
      DataStore::self_ = new InMemoryDataStore(0, next);
    }
  }

  void DataStore::finalize() {
    if (DataStore::self_ != NULL) {
      delete DataStore::self_;
      DataStore::self_ = NULL;
    }
  }

  void DataStore::zeroize() {
    Zeroizer zero;
    vector<long> dlst;
    for (size_t i = 0; i < value_[0]; i++) {
      dlst.push_back(i);
    }
    load_integers(dlst, zero);
  }

  ////////////////////////////////////////////////////////////////////////////
  // Method implementation for class DataStoreCommonImpl
  ////////////////////////////////////////////////////////////////////////////

  void DataStoreCommonImpl::perform_duplicate(DataStore* ds) {
    class Copier : public Mapper {
    public:
      int operator()(DataStore* inds, DataStore* outds,
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
    ds->set(value_);
    View view(size_);
    for (size_t i = 0; i < size_; i++) {
      view.set_dim(i, View::SplitNone);
    }
    map(copier, view, ds);
  }

  size_t DataStoreCommonImpl::key_to_index(const Key& key) {
    size_t idx = 0;
    for (size_t i = 0; i < size_; i++) {
      idx += key.dim(i) * icache_.dim_offset(i);
    }
    return idx;
  }

  Key DataStoreCommonImpl::index_to_key(const size_t index) {
    return icache_.i2k(index);
  }

  size_t DataStoreCommonImpl::key_to_viewed_index(const Key& key,
						  const View& view) {
    // Currently, It returns column-ordered index
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

  Key DataStoreCommonImpl::key_to_viewed_key(const Key& key, const View& view)
  {
    size_t viewed_key_size = 0;
    for (size_t i = 0; i < size_; i++) {
      if (view.dim(i) != View::SplitNone) {
	viewed_key_size += 1;
      }
    }

    Key viewed_key(viewed_key_size);
    size_t viewed_key_idx = 0;
    for (size_t i = 0; i < size_; i++) {
      if (view.dim(i) == View::SplitAll) {
	viewed_key.set_dim(viewed_key_idx, key.dim(i));
	viewed_key_idx += 1;
      } else if (view.dim(i) > 0) {  // split size is specified
	size_t blk_siz = value_[i] / view.dim(i);
	viewed_key.set_dim(viewed_key_idx, key.dim(i) / blk_siz);
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
      if (view.dim(i) == View::SplitAll) {
	viewed_key.set_dim(i, key.dim(i));
      } else if (view.dim(i) > 0) {  // split size is specified
	size_t blk_siz = value_[i] / view.dim(i);
	viewed_key.set_dim(i, key.dim(i) / blk_siz);
      } else {
	viewed_key.set_dim(i, 0);
      }
    }
    return viewed_key;
  }
#endif

  void DataStoreCommonImpl::check_view(const View& view) {
    if (size_ != view.size()) {
      throw runtime_error("Dimension size of the DataStore and view "
			  "should be same.");
    }
    for (size_t i = 0; i < size_; i++) {
      if (!(view.dim(i) == View::SplitAll || view.dim(i) == View::SplitNone)) {
	if (value_[i] % view.dim(i) != 0) {
	  ostringstream os;
	  os << "View/Split could not divide the DataStore." << to_string();
	  throw runtime_error(os.str());
	}
      }
    }
  }

  void DataStoreCommonImpl::check_key_range(const Key& key) {
    if (size_ != key.size()) {
      throw runtime_error("Dimension size of Key should be same as "
			  "that of DataStore.");
    }
    for (size_t i = 0; i < size_; i++) {
      if (key.dim(i) >= value_[i]) {
	ostringstream os;
	os << "Dimension " << (i+1) << " of Key" << key.to_string()
	   << " is out of range.";
	throw runtime_error(os.str());
      }
    }
  }

  void DataStoreCommonImpl::check_map_args(const View& view, DataStore* outds)
  {
    if (outds == NULL) {
      throw runtime_error("The output DataStore should not be NULL.");
    }
    if (size_ != view.size()) {
      throw runtime_error("Dimension size of the input DataStore and "
			  "view should be same.");
    }
    for (size_t i = 0; i < size_; i++) {
      if (!(view.dim(i) == View::SplitAll || view.dim(i) == View::SplitNone)) {
	if (value_[i] % view.dim(i) != 0) {
	  ostringstream os;
	  os << "View could not divide the input DataStore." << to_string();
	  throw runtime_error(os.str());
	}
      }
    }
  }

  ////////////////////////////////////////////////////////////////////////////
  // Method implementation for class IndexCache
  ////////////////////////////////////////////////////////////////////////////

  IndexCache::IndexCache()
    : i2k_table_(vector<Key>()), doffset_table_(vector<size_t>()) {}

  void IndexCache::initialize(const size_t* sizes, const size_t i2k_len,
			      const size_t dim_siz) {
    doffset_table_.reserve(dim_siz);
    for (size_t i = 0; i < dim_siz; i++) {
      doffset_table_[i] = 1;
      for (size_t j = i+1; j < dim_siz; j++) {
	doffset_table_[i] *= sizes[j];
      }
    }

    i2k_table_.reserve(i2k_len);
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (size_t index = 0; index < i2k_len; index++) {
      Key key(dim_siz);
      size_t _index = index;
      for (size_t i = 0; i < dim_siz; i++) {
	key.set_dim(i, _index / doffset_table_[i]);
	_index %= doffset_table_[i];
      }
      i2k_table_[index] = key;
    }
  }

  Key IndexCache::i2k(const size_t index) const {
    return i2k_table_[index];
  }

  size_t IndexCache::dim_offset(const size_t dim) const {
    return doffset_table_[dim];
  }

  ////////////////////////////////////////////////////////////////////////////
  // Method implementation for class SimpleFileDataStore
  ////////////////////////////////////////////////////////////////////////////

  SimpleFileDataStore::~SimpleFileDataStore() {
    string fname = filename();
    delete_file(fname);
  }

  void SimpleFileDataStore::add(const Key& key, const Data& data) {
    base::add(key, data);
    data_updated_ = true;
  }

  void SimpleFileDataStore::set_from(const vector<DataStore*>& dslist) {
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

  void SimpleFileDataStore::split_to(vector<DataStore*>& dslist) {
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

  void SimpleFileDataStore::map(Mapper& m, const View& view, DataStore* outds)
  {
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

  DataStore* SimpleFileDataStore::duplicate() {
    SimpleFileDataStore* ds = new SimpleFileDataStore(size_, kmrnext_);
    perform_duplicate(ds);
    return ds;
  }

  bool SimpleFileDataStore::store() {
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

  bool SimpleFileDataStore::load() {
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

  void SimpleFileDataStore::clear_cache() {
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

}

#ifdef BACKEND_SERIAL
#include "data_store_serial.cpp"
#elif defined BACKEND_KMR
#include "data_store_kmr.cpp"
#endif
