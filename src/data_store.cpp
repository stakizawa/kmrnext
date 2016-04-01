#include "../config.hpp"
#include "kmrnext.hpp"

namespace kmrnext {

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
