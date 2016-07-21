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
    }
  }

  KMRNext::IOMode DataStore::io_mode() {
    if (kmrnext_ == NULL) {
      throw runtime_error("KMRNext instance should be set to a DataStore.");
    }
    return kmrnext_->io_mode();
  }

  void DataStore::load_files(const vector<string>& files,
			     Loader<string>& loader) {
    load_array(files, loader, kmrnext_, this, value_, size_);
  }

  void DataStore::load_integers(const vector<long>& ints,
				Loader<long>& loader) {
    load_array(ints, loader, kmrnext_, this, value_, size_);
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

  bool DataStore::store() {
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
    size_t cur_buf_siz = kDefaultWriteBufferSize;
    char *buf = static_cast<char*>(calloc(cur_buf_siz, sizeof(char)));

    for (size_t i = 0; i < dlist_size_; i++) {
      if (!dlist_[i].is_set()) {
	continue;
      }
      Data *d = dlist_[i].data();
      size_t d_siz = d->size();
      if (d_siz == 0) {
	continue;
      }
      void  *d_val = d->value();
      size_t buf_siz = d_siz;
      if (buf_siz > cur_buf_siz) {
	cur_buf_siz = buf_siz;
	buf = static_cast<char*>(realloc(buf, cur_buf_siz));
      }
      // TODO maybe memcpy() is not necessary
      memcpy(buf, d_val, d_siz);
      fout.write(buf, static_cast<streamsize>(buf_siz));
      dlist_[i].written(write_offset, buf_siz);
      write_offset += buf_siz;
    }
    free(buf);
    fout.flush();
    fout.close();
    data_updated_ = false;
    data_cached_ = true;
    return true;
  }

  bool DataStore::load() {
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
      dlist_[i].restore(buf);
    }
    free(buf);
    data_cached_ = true;
    return true;
  }

  void DataStore::clear_cache() {
    if (!data_cached_) {
      return;
    }
#ifdef _OPENMP
    #pragma omp parallel for
#endif
    for (size_t i = 0; i < dlist_size_; i++) {
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
    split_ds0.set_dim(0, true);
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
    v.set_dim(0, true);
    ds0->map(wloader, v, ds);
    delete ds0;

#ifdef BACKEND_KMR
    View split_ds(ds_dims_siz);
    if (array.size() == 1) {
      // Use Split <T, F, ..> so that data will be distributed to nodes
      // whose count is the top most dimension of the DataStore.
      split_ds.set_dim(0, true);
      for (size_t i = 1; i < ds_dims_siz; i++) {
	split_ds.set_dim(i, false);
      }
    } else {
      // Use Split <T, .., T, F, ..> so that each node stores
      // the contents of each array element.
      bool vval = true;
      size_t remain = array.size();
      for (size_t i = 0; i < ds_dims_siz; i++) {
	split_ds.set_dim(i, vval);
	remain /= ds->dim(i);
	if (remain == 1) {
	  vval = false;
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
