#include "../config.hpp"
#include "kmrnext.hpp"

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include "util.hpp"

namespace {
  using namespace std;
  using namespace kmrnext;

  const size_t kDefaultWriteBufferSize = 1048576;

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

  DataStore *DataStore::self_;

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

  void DataStore::set(const size_t *val) {
    if (dlist_size_ != 0) {
      throw runtime_error("DataStore is already initialized.");
    }

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
      dlist_[i] = new DataElement();
    }
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
	    throw runtime_error("Data count of one of DataStore is wrong.");
	  }
	}
      }
    }

    size_t sizes[kMaxDimensionSize];
    sizes[0] = dslist.size();
    DataStore *ds0 = dslist.at(0);
    for (size_t i = 1; i < size_; i++) {
      sizes[i] = ds0->value_[i-1];
    }
    set(sizes);

    size_t offset = 0;
    for (size_t i = 0; i < dslist.size(); i++) {
      DataStore *src = dslist.at(i);
      for (size_t j = 0; j < src->dlist_size_; j++) {
	dlist_[offset + j]->set(src->dlist_[j]->data());
      }
      offset += src->dlist_size_;
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

    size_t split_dims[kMaxDimensionSize];
    for (size_t i = 1; i < size_; i++) {
      split_dims[i-1] = value_[i];
    }

    size_t offset = 0;
    for (size_t i = 0; i < dslist.size(); i++) {
      DataStore *dst = dslist.at(i);
      dst->set(split_dims);
      for (size_t j = 0; j < dst->dlist_size_; j++) {
	dst->dlist_[j]->set(dlist_[offset + j]->data());
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
    // Currently, It returns column-ordered index
#if 1
    // It returns column-ordered index
    size_t idx = 0;
    for (long i = static_cast<long>(size_) - 1; i >= 0; i--) {
      if (view.dim(i) == View::SplitAll) {
	size_t offset = 1;
	for (long j = i-1; j >= 0; j--) {
	  if (view.dim(j) == View::SplitAll) {
	    offset *= value_[j];
	  }
	}
	idx += key.dim(i) * offset;
      }
    }
    return idx;
#else
    // It returns row-ordered index
    size_t idx = 0;
    for (size_t i = 0; i < size_; i++) {
      if (view.dim(i) == View::SplitAll) {
	size_t offset = 1;
	for (size_t j = i+1; j < size_; j++) {
	  if (view.dim(j) == View::SplitAll) {
	    offset *= value_[j];
	  }
	}
	idx += key.dim(i) * offset;
      }
    }
    return idx;
#endif
  }

  Key DataStore::key_to_viewed_key(const Key& key, const View& view) {
    size_t viewed_key_size = 0;
    for (size_t i = 0; i < size_; i++) {
      if (view.dim(i) == View::SplitAll) {
	viewed_key_size += 1;
      }
    }

    Key viewed_key(viewed_key_size);
    size_t viewed_key_idx = 0;
    for (size_t i = 0; i < size_; i++) {
      if (view.dim(i) == View::SplitAll) {
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
      if (view.dim(i) == View::SplitAll) {
	viewed_key.set_dim(i, key.dim(i));
      } else {
	viewed_key.set_dim(i, 0);
      }
    }
    return viewed_key;
  }
#endif

  DataElement* DataStore::data_element_at(const Key& key) {
    size_t idx = key_to_index(key);
    return dlist_[idx];
  }

  void DataStore::check_view(const View& view) {
    if (size_ != view.size()) {
      throw runtime_error("Dimension size of the DataStore and view "
			  "should be same.");
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
  }

  SimpleFileDataStore::~SimpleFileDataStore() {
    string fname = filename();
    delete_file(fname);
  }

  void SimpleFileDataStore::set(const size_t *val) {
    if (dlist_size_ != 0) {
      throw runtime_error("DataStore is already initialized.");
    }

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
      dlist_[i] = new SimpleFileDataElement();
    }
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
      SimpleFileDataStore *ds =
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
      if (!dlist_[i]->is_set()) {
	continue;
      }
      Data *d = dlist_[i]->data();
      size_t d_siz = d->size();
      if (d_siz == 0) {
	continue;
      }
      char *d_val = static_cast<char*>(d->value());
      fout.write(d_val, static_cast<streamsize>(d_siz));
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
    char *buf = static_cast<char*>(calloc(file_siz, sizeof(char)));
    ifstream fin;
    fin.open(fname.c_str(), ios::in|ios::binary);
    fin.read(buf, static_cast<streamsize>(file_siz));
    fin.close();

#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (size_t i = 0; i < dlist_size_; i++) {
      dynamic_cast<SimpleFileDataElement*>(dlist_[i])->restore(buf);
    }
    free(buf);
    data_cached_ = true;
    return true;
  }

  void SimpleFileDataStore::clear_cache() {
    if (!data_cached_) {
      return;
    }
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (size_t i = 0; i < dlist_size_; i++) {
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

    DataStore* ds0 = new DataStore(1, next);
    ds0->set_dim(0, array.size());
    Key key(1);
    for (size_t i = 0; i < array.size(); i++) {
      key.set_dim(0, i);
      char *buf;
      size_t buf_siz;
      serialize(array.at(i), &buf, &buf_siz);
      Data dat(buf, buf_siz);
      ds0->add(key, dat);
    }
#ifdef BACKEND_KMR
    // Use Split <T> so that each process loads an array element.
    View split_ds0(1);
    split_ds0.set_dim(0, View::SplitAll);
    ds0->set_split(split_ds0);
#endif

    // Define a mapper for the loader
    class WrappedLoader : public DataStore::Mapper {
    public:
      DataStore::Loader<T>& loader_;

      WrappedLoader(DataStore::Loader<T>& ldr) : loader_(ldr) {}
      int operator()(DataStore *inds, DataStore *outds,
		     Key& k, vector<DataPack>& dps,
		     DataStore::MapEnvironment& env)
      {
	T *val;
	deserialize(static_cast<char*>(dps.at(0).data().value()),
		    dps.at(0).data().size(), &val);
	loader_(outds, *val);
	delete val;
	return 0;
      }
    } wloader(loader);

    View v(1);
    v.set_dim(0, View::SplitAll);
    ds0->map(wloader, v, ds);
    delete ds0;

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
