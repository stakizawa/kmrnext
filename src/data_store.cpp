#include "../config.hpp"
#include "kmrnext.hpp"

#include <iostream>
#include <fstream>
#include <cstring>
#include "util.hpp"

#ifdef _OPENMP
#define OMP_FOR_CHUNK_SIZE 4
#endif

namespace {
  using namespace std;
  using namespace kmrnext;

  // A data loader class that sets long integer 0 to all data elements.
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

  template <typename T>
  void load_array(const vector<T>& array, DataStore::Loader<T>& loader,
		  KMRNext* next, DataStore* ds,
		  size_t* ds_dims, size_t ds_dims_siz);

  // It serializes a string.
  void serialize(const string& str, char** buf, size_t* buf_siz);

  // It deserializes a string.
  void deserialize(char* buf, size_t buf_siz, string** str);

  // It serializes an integer.
  void serialize(const long& val, char** buf, size_t* buf_siz);

  // It deserializes an integer.
  void deserialize(char* buf, size_t buf_siz, long** val);
}

namespace kmrnext {

  DataStore* DataStore::self_;

  void DataStore::initialize(KMRNext* next) {
    if (DataStore::self_ == NULL) {
      DataStore::self_ = new DataStore(0, next);
    }
  }

  void DataStore::finalize() {
    if (DataStore::self_ != NULL) {
      delete DataStore::self_;
      DataStore::self_ = NULL;
    }
  }

  void DataStore::set(const size_t* val) {
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
    dlist_.reserve(dlist_size_);
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (size_t i = 0; i < dlist_size_; i++) {
      dlist_[i] = NULL;
    }

    icache_.initialize(value_, dlist_size_, size_);
  }

  void DataStore::set_dim(const size_t idx, const size_t siz) {
#if VALIDATION
    if (dlist_size_ != 0) {
      throw runtime_error("DataStore is already initialized.");
    }
#endif
    base::set_dim(idx, siz);
  }

  void DataStore::zeroize() {
    Zeroizer zero;
    vector<long> dlst;
    for (size_t i = 0; i < value_[0]; i++) {
      dlst.push_back(i);
    }
    load_integers(dlst, zero);
  }

  void DataStore::set_from(const vector<DataStore*>& dslist) {
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
	DataStore* src = dslist.at(i);
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
    DataStore* ds0 = dslist.at(0);
    for (size_t i = 1; i < size_; i++) {
      sizes[i] = ds0->value_[i-1];
    }
    set(sizes);

    size_t offset = 0;
    for (size_t i = 0; i < dslist.size(); i++) {
      DataStore* src = dslist.at(i);
      for (size_t j = 0; j < src->dlist_size_; j++) {
	if (dlist_[offset + j] == NULL) {
	  dlist_[offset + j] = __create_de();
	}
	if (src->dlist_[j] == NULL) {
	  continue;
	}
	dlist_[offset + j]->set(src->dlist_[j]->value(),
				src->dlist_[j]->size());
      }
      offset += src->dlist_size_;
    }
  }

  void DataStore::split_to(vector<DataStore*>& dslist) {
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
	if (expected_dim_size != (*itr)->size_) {
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
      DataStore* dst = dslist.at(i);
      dst->set(split_dims);
      for (size_t j = 0; j < dst->dlist_size_; j++) {
	if (dst->dlist_[j] == NULL) {
	  dst->dlist_[j] = __create_de();
	}
	if (dlist_[offset + j] == NULL) {
	  continue;
	}
	dst->dlist_[j]->set(dlist_[offset + j]->value(),
			    dlist_[offset + j]->size());
      }
      offset += dst->dlist_size_;
    }
  }

  DataStore* DataStore::duplicate() {
    DataStore* ds = new DataStore(size_, kmrnext_);
    __duplicate(ds);
    return ds;
  }

  void DataStore::__duplicate(DataStore* ds) {
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

  void DataStore::load_files(const vector<string>& files,
			     Loader<string>& loader) {
    load_array(files, loader, kmrnext_, this, value_, size_);
  }

  void DataStore::load_integers(const vector<long>& ints,
				Loader<long>& loader) {
    load_array(ints, loader, kmrnext_, this, value_, size_);
  }

  DataElement* DataStore::__create_de() {
    return new DataElement();
  }

  size_t DataStore::key_to_index(const Key& key) {
    size_t idx = 0;
    for (size_t i = 0; i < size_; i++) {
      idx += key.dim(i) * icache_.dim_offset(i);
    }
    return idx;
  }

  Key DataStore::index_to_key(const size_t index) {
    return icache_.i2k(index);
  }

  size_t DataStore::key_to_viewed_index(const Key& key, const View& view) {
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

  Key DataStore::key_to_viewed_key(const Key& key, const View& view) {
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

  DataElement* DataStore::data_element_at(const Key& key) {
    size_t idx = key_to_index(key);
    if (dlist_[idx] == NULL) {
      dlist_[idx] = __create_de();
    }
    return dlist_[idx];
  }

  void DataStore::check_view(const View& view) {
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

  void DataStore::check_key_range(const Key& key) {
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

  void DataStore::check_map_args(const View& view, DataStore* outds) {
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

  DataStore::IndexCache::IndexCache() : i2k_table_(vector<Key>()),
					doffset_table_(vector<size_t>()) {}

  void DataStore::IndexCache::initialize(const size_t* sizes,
					 const size_t i2k_len,
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

  Key DataStore::IndexCache::i2k(const size_t index) const {
    return i2k_table_[index];
  }

  size_t DataStore::IndexCache::dim_offset(const size_t dim) const {
    return doffset_table_[dim];
  }

  SimpleFileDataStore::~SimpleFileDataStore() {
    string fname = filename();
    delete_file(fname);
  }

  void SimpleFileDataStore::add(const Key& key, const Data& data) {
    DataStore::add(key, data);
    data_updated_ = true;
  }

  void SimpleFileDataStore::set_from(const vector<DataStore*>& dslist) {
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (size_t i = 0; i < dslist.size(); i++) {
      dynamic_cast<SimpleFileDataStore*>(dslist.at(i))->load();
    }
    DataStore::set_from(dslist);
    store();
    clear_cache();
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (size_t i = 0; i < dslist.size(); i++) {
      dynamic_cast<SimpleFileDataStore*>(dslist.at(i))->clear_cache();
    }
  }

  void SimpleFileDataStore::split_to(vector<DataStore*>& dslist) {
    load();
    DataStore::split_to(dslist);
    for (size_t i = 0; i < dslist.size(); i++) {
      SimpleFileDataStore* ds =
	dynamic_cast<SimpleFileDataStore*>(dslist.at(i));
      ds->store();
      ds->clear_cache();
    }
    clear_cache();
  }

  void SimpleFileDataStore::map(Mapper& m, const View& view, DataStore* outds)
  {
    store();
    load();
    DataStore::map(m, view, outds);
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
    __duplicate(ds);
    return ds;
  }

  DataElement* SimpleFileDataStore::__create_de() {
    return new SimpleFileDataElement();
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
      if (dlist_[i] == NULL || !dlist_[i]->is_set()) {
	continue;
      }
      size_t d_siz = dlist_[i]->size();
      if (d_siz == 0) {
	continue;
      }
      fout.write(static_cast<char*>(dlist_[i]->value()),
		 static_cast<streamsize>(d_siz));
      dynamic_cast<SimpleFileDataElement*>(dlist_[i])->written(write_offset,
							       d_siz);
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
      if (dlist_[i] == NULL) {
	continue;
      }
      dynamic_cast<SimpleFileDataElement*>(dlist_[i])->restore(buf);
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
      if (dlist_[i] == NULL) {
	continue;
      }
      dynamic_cast<SimpleFileDataElement*>(dlist_[i])->clear_cache();
    }
    data_cached_ = false;
  }

}

#ifdef BACKEND_SERIAL
#include "data_store_serial.cpp"
#elif defined BACKEND_KMR
#include "data_store_kmr.cpp"
#endif

namespace {

  template <typename T>
  void load_array(const vector<T>& array, DataStore::Loader<T>& loader,
		  KMRNext* next, DataStore* ds,
		  size_t* ds_dims, size_t ds_dims_siz) {
#if VALIDATION
    // Check if the size of array is same as the multiple of dimension.
    {
      size_t prod = 1;
      for (size_t i = 0; i < ds_dims_siz; i++) {
	prod *= ds_dims[i];
	if (prod == array.size()) {
	  break;
	}
      }
      if (array.size() != 1 && array.size() != prod) {
	throw runtime_error("The size of array should be 1 or match the "
			    "product of dimension sizes of the DataStore.");
      }
    }
#endif

    DataStore ds0(1, next);
    ds0.set_dim(0, array.size());
    Key key(1);
    for (size_t i = 0; i < array.size(); i++) {
      key.set_dim(0, i);
      char* buf;
      size_t buf_siz;
      serialize(array.at(i), &buf, &buf_siz);
      Data dat(buf, buf_siz);
      ds0.add(key, dat);
    }
#ifdef BACKEND_KMR
    // Use Split <T> so that each process loads an array element.
    View split_ds0(1);
    split_ds0.set_dim(0, View::SplitAll);
    ds0.set_split(split_ds0);
#endif

    // Define a mapper for the loader
    class WrappedLoader : public DataStore::Mapper {
    public:
      DataStore::Loader<T>& loader_;

      WrappedLoader(DataStore::Loader<T>& ldr) : loader_(ldr) {}
      int operator()(DataStore* inds, DataStore* outds,
		     Key& k, vector<DataPack>& dps,
		     DataStore::MapEnvironment& env)
      {
	T* val;
	deserialize(static_cast<char*>(dps.at(0).data().value()),
		    dps.at(0).data().size(), &val);
	loader_(outds, *val);
	delete val;
	return 0;
      }
    } wloader(loader);

    View v(1);
    v.set_dim(0, View::SplitAll);
    ds0.map(wloader, v, ds);

#ifdef BACKEND_KMR
    View split_ds(ds_dims_siz);
    if (array.size() == 1) {
      // Use Split <All, None, ..> so that data will be distributed to nodes
      // whose count is the top most dimension of the DataStore.
      split_ds.set_dim(0, View::SplitAll);
      for (size_t i = 1; i < ds_dims_siz; i++) {
	split_ds.set_dim(i, View::SplitNone);
      }
    } else {
      // Use Split <All, .., All, None, ..> so that each node stores
      // the contents of each array element.
      long vval = View::SplitAll;
      size_t remain = array.size();
      for (size_t i = 0; i < ds_dims_siz; i++) {
	split_ds.set_dim(i, vval);
	remain /= ds->dim(i);
	if (remain == 1) {
	  vval = View::SplitNone;
	}
      }
    }
    ds->set_split(split_ds);
#endif
  }

  void serialize(const string& str, char** buf, size_t* buf_siz) {
    *buf = const_cast<char*>(str.c_str());
    *buf_siz = str.size() + 1; // +1 for '\0'
  }

  void deserialize(char* buf, size_t buf_siz, string** str) {
    *str = new string(buf);
  }

  void serialize(const long& val, char** buf, size_t* buf_siz) {
    long* v = const_cast<long*>(&val);
    *buf = reinterpret_cast<char*>(v);
    *buf_siz = sizeof(long);
  }

  void deserialize(char* buf, size_t buf_siz, long** val) {
    *val = new long;
    **val = *reinterpret_cast<long*>(buf);
  }

}
