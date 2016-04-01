#include "../config.hpp"
#include "kmrnext.hpp"

namespace kmrnext {

  void DataStore::load_files(const vector<string>& files,
			     Loader<string>& loader) {
    // Check if the size of array is same as the multiple of dimension.
    {
      size_t prod = 1;
      for (size_t i = 0; i < size_; i++) {
	prod *= value_[i];
	if (prod == files.size()) {
	  break;
	}
      }
      if (files.size() != 1 && files.size() != prod) {
	throw runtime_error("The size of array should be 1 or match the "
			    "product of dimension sizes of the DataStore.");
      }
    }

    // Create a filename DataStore
    DataStore* ds0 = new DataStore(1, kmrnext_);
    ds0->set_dim(0, files.size());
    Key key(1);
    for (size_t i = 0; i < files.size(); i++) {
      key.set_dim(0, i);
      void *cstr = (void*)files.at(i).c_str();
      size_t dat_siz = files.at(i).size() + 1;
      Data dat(cstr, dat_siz);
      ds0->add(key, dat);
    }
#ifdef BACKEND_KMR
    // TODO set physical_view to ds0
    throw runtime_error("Not implemented yet.");
#endif

    // Define a mapper for the loader
    class WrappedLoader : public Mapper {
    public:
      Loader<string>& loader_;

      WrappedLoader(Loader<string>& loader) : loader_(loader) {}
      int operator()(DataStore *inds, DataStore *outds,
		     Key& key, vector<DataPack>& dps,
		     MapEnvironment& env)
      {
	char *val = (char*)dps.at(0).data()->value();
	string filename(val);
	loader_(outds, filename);
	return 0;
      }
    } wloader(loader);

    View v(1);
    v.set_dim(0, true);
    ds0->map(wloader, v, this);
    delete ds0;

#ifdef BACKEND_KMR
    // TODO set physical_view to this
    throw runtime_error("Not implemented yet.");
#endif
  }

  void DataStore::load_integers(const vector<int>& ints,
				Loader<int>& loader) {
    // Check if the size of array is same as the multiple of dimension.
    {
      size_t prod = 1;
      for (size_t i = 0; i < size_; i++) {
	prod *= value_[i];
	if (prod == ints.size()) {
	  break;
	}
      }
      if (ints.size() != 1 && ints.size() != prod) {
	throw runtime_error("The size of array should be 1 or match the "
			    "product of dimension sizes of the DataStore.");
      }
    }

    // Create an integer DataStore
    DataStore* ds0 = new DataStore(1, kmrnext_);
    ds0->set_dim(0, ints.size());
    Key key(1);
    for (size_t i = 0; i < ints.size(); i++) {
      key.set_dim(0, i);
      Data dat((void*)&(ints.at(i)), sizeof(int));
      ds0->add(key, dat);
    }
#ifdef BACKEND_KMR
    // TODO set physical_view to ds0
    throw runtime_error("Not implemented yet.");
#endif

    // Define a mapper for the loader
    class WrappedLoader : public Mapper {
    public:
      Loader<int>& loader_;

      WrappedLoader(Loader<int>& loader) : loader_(loader) {}
      int operator()(DataStore *inds, DataStore *outds,
		     Key& key, vector<DataPack>& dps,
		     MapEnvironment& env)
      {
	int val = *(int*)dps.at(0).data()->value();
	loader_(outds, val);
	return 0;
      }
    } wloader(loader);

    View v(1);
    v.set_dim(0, true);
    ds0->map(wloader, v, this);
    delete ds0;

#ifdef BACKEND_KMR
    // TODO set physical_view to this
    throw runtime_error("Not implemented yet.");
#endif
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
