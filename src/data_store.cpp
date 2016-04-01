#include "../config.hpp"
#include "kmrnext.hpp"

#include <iostream>

namespace {
  using namespace std;
  using namespace kmrnext;

  template <typename T>
  void load_array(const vector<T>& array, DataStore::Loader<T>& loader,
		  KMRNext* next, DataStore* ds,
		  size_t* ds_dims, size_t ds_dims_siz);


  // It serializes a string.
  void serialize(const string& str, char** buf, size_t* buf_siz);

  // It deserializes a string.
  void deserialize(char* buf, size_t buf_siz, string** str);

  // It serializes an integer.
  void serialize(const int& val, char** buf, size_t* buf_siz);

  // It deserializes an integer.
  void deserialize(char* buf, size_t buf_siz, int** val);
}

namespace kmrnext {

  void DataStore::load_files(const vector<string>& files,
			     Loader<string>& loader) {
    load_array(files, loader, kmrnext_, this, value_, size_);
  }

  void DataStore::load_integers(const vector<int>& ints,
				Loader<int>& loader) {
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

  void DataStore::check_collate_args(const View& view) {
    if (size_ != view.size()) {
      throw runtime_error("Dimension size of the DataStore and "
			  "view should be same.");
    }
    if (data_size_ == 0) {
      throw runtime_error("Data should be set to the DataStore.");
    }
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
    // TODO set physical_view to ds0
    throw runtime_error("Not implemented yet.");
#endif

    // Define a mapper for the loader
    class WrappedLoader : public DataStore::Mapper {
    public:
      DataStore::Loader<T>& loader_;

      WrappedLoader(DataStore::Loader<T>& loader) : loader_(loader) {}
      int operator()(DataStore *inds, DataStore *outds,
		     Key& key, vector<DataPack>& dps,
		     DataStore::MapEnvironment& env)
      {
	T *val;
	deserialize((char*)dps.at(0).data()->value(),
		    dps.at(0).data()->size(), &val);
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
    // TODO set physical_view to this
    throw runtime_error("Not implemented yet.");
#endif
  }

  void serialize(const string& str, char** buf, size_t* buf_siz) {
    *buf = (char*)str.c_str();
    *buf_siz = str.size() + 1; // +1 for '\0'
  }

  void deserialize(char* buf, size_t buf_siz, string** str) {
    *str = new string(buf);
  }

  void serialize(const int& val, char** buf, size_t* buf_siz) {
    *buf = (char*)&val;
    *buf_siz = sizeof(int);
  }

  void deserialize(char* buf, size_t buf_siz, int** val) {
    *val = new int[1];
    **val = (int)*buf;
  }

}
