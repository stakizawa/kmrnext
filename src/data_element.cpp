#include "../config.hpp"

#include "kmrnext.hpp"

namespace kmrnext {

#ifdef BACKEND_SERIAL
  DataElement::DataElement()
    : data_(NULL), data_set_(false),
      data_updated_(false), data_file_offset_(0), data_file_size_(0) {}
#elif BACKEND_KMR
  DataElement::DataElement()
    : data_(NULL), data_set_(false),
      data_updated_(false), data_file_offset_(0), data_file_size_(0),
      owner_(-1), shared_(false) {}
#endif

  DataElement::~DataElement() {
    if (data_set_) {
      delete data_;
    }
  }

  void DataElement::set(const Data* dat) {
    set_data(dat);
  }

  void DataElement::replace(const Data* dat) {
    set_data(dat, true);
  }

  void DataElement::set_data(const Data* dat, bool overwrite) {
    if (dat == NULL) {
      return;
    }
    if (overwrite) {
      if (data_set_) {
	delete data_;
      }
    } else {
      if (data_set_) {
	throw runtime_error("Data is already set value.");
      }
    }
    data_ = new Data(dat->value(), dat->size());
    data_->allocate();
    data_set_ = true;

    data_updated_ = true;
  }

  void DataElement::clear() {
    if (data_set_) {
      delete data_;
    }
    data_ = NULL;
    data_set_ = false;

    data_updated_ = false;
    data_file_offset_ = 0;
    data_file_size_ = 0;

#ifdef BACKEND_KMR
    owner_ = -1;
    shared_ = false;
#endif
  }

  void DataElement::restore(char *buf) {
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
    if (data_ != NULL) {
      throw runtime_error("Data should be NULL.");
    }
    data_ = new Data(buf + data_file_offset_, data_file_size_);
    data_->allocate();
  }

  void DataElement::written(size_t start_pos, size_t written_siz) {
    data_updated_ = false;
    data_file_offset_ = start_pos;
    data_file_size_ = written_siz;
  }

  void DataElement::clear_cache() {
    if (data_ != NULL) {
      delete data_;
      data_ = NULL;
    }
    data_updated_ = false;
  }

}
