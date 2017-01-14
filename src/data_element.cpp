#include "../config.hpp"

#include <cstdlib>
#include <cstring>
#include "kmrnext.hpp"

namespace kmrnext {

#ifdef BACKEND_SERIAL
  DataElement::DataElement()
    : value_(NULL), value_size_(0), data_set_(false) {}
#elif BACKEND_KMR
  DataElement::DataElement()
    : value_(NULL), value_size_(0), data_set_(false),
      owner_(-1), shared_(false) {}
#endif

  DataElement::~DataElement() {
    if (data_set_) {
      free(value_); // TODO delete[]?
    }
  }

  void DataElement::set(const void* data_value, const size_t data_size) {
    set_data(data_value, data_size);
  }

  void DataElement::replace(const void* data_value, const size_t data_size) {
    set_data(data_value, data_size, true);
  }

  void DataElement::set_data(const void* val, const size_t siz,
			     bool overwrite) {
    if (val == NULL) {
      return;
    }
    if (overwrite) {
      if (data_set_) {
	free(value_);  // TODO delete[]?
      }
    } else {
#if VALIDATION
      if (data_set_) {
	throw runtime_error("Data is already set value.");
      }
#endif
    }
    value_size_ = siz;
    value_ = static_cast<void*>(calloc(value_size_, sizeof(char)));
    memcpy(value_, val, value_size_);
    data_set_ = true;
  }

  void DataElement::clear() {
    if (data_set_) {
      value_size_ = 0;
      free(value_);  // TODO delete[]?
      value_ = NULL;
    }
    data_set_ = false;
#ifdef BACKEND_KMR
    owner_ = -1;
    shared_ = false;
#endif
  }

  SimpleFileDataElement::SimpleFileDataElement()
    : base(),
      data_updated_(false), data_file_offset_(0), data_file_size_(0) {}

  void SimpleFileDataElement::set_data(const void* val, const size_t siz,
				       bool overwrite) {
    DataElement::set_data(val, siz, overwrite);
    data_updated_ = true;
  }

  void SimpleFileDataElement::clear() {
    DataElement::clear();
    data_updated_ = false;
    data_file_offset_ = 0;
    data_file_size_ = 0;
  }

  void SimpleFileDataElement::restore(char* buf) {
    if (!data_set_) {
      // Skip as Data is removed
      return;
    }
    if (data_updated_) {
      // Skip as Data is updated in memory
      return;
    }
    if (data_file_size_ == 0) {
      throw runtime_error("Data can not be found in the file.");
    }
    if (value_ != NULL) {
      throw runtime_error("Data should not be set.");
    }
    value_size_ = data_file_size_;
    value_ = static_cast<void*>(calloc(value_size_, sizeof(char)));
    memcpy(value_, buf + data_file_offset_, value_size_);
  }

  void SimpleFileDataElement::written(size_t start_pos, size_t written_siz) {
    data_updated_ = false;
    data_file_offset_ = start_pos;
    data_file_size_ = written_siz;
  }

  void SimpleFileDataElement::clear_cache() {
    if (value_ != NULL) {
      value_size_ = 0;
      free(value_);  // TODO delete[]?
      value_ = NULL;
    }
    data_updated_ = false;
  }

}
