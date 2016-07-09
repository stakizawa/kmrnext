#include "../config.hpp"

#include <cstdlib>
#include <cstring>
#include "kmrnext.hpp"

namespace kmrnext {

#ifdef BACKEND_SERIAL
  Data::Data() :
    value_(NULL), value_size_(0), value_allocated_(false),
    data_set_(false), data_file_offset_(0), data_file_size_(0) {}

  Data::Data(void *val, const size_t val_siz) :
    value_(val), value_size_(val_siz), value_allocated_(false),
    data_set_(true), data_file_offset_(0), data_file_size_(0) {}
#elif defined BACKEND_KMR
  Data::Data() :
    value_(NULL), value_size_(0), value_allocated_(false),
    data_set_(false), data_file_offset_(0), data_file_size_(0),
    owner_(-1), shared_(false) {}

  Data::Data(void *val, const size_t val_siz) :
    value_(val), value_size_(val_siz), value_allocated_(false),
    data_set_(true), data_file_offset_(0), data_file_size_(0),
    owner_(-1), shared_(false) {}
#endif

  Data::~Data() {
    if (value_allocated_) {
      free(value_);
    }
  }

  void Data::set_value(const Data& src) {
    if (data_set_) {
      throw runtime_error("Data is already set value.");
    }
    copy_deep(src);
  }

  void Data::update_value(const Data& src) {
    if (!data_set_) {
      throw runtime_error("Data is not set value.");
    }
    copy_deep(src, true);
  }

  void Data::copy_deep(const Data& src, bool overwrite) {
    if (overwrite) {
      if (src.value_ == NULL) {
	throw runtime_error("The copy target Data should not be NULL.");
      }
      if (value_ == NULL) {
	value_ = static_cast<void*>(calloc(src.value_size_, sizeof(char)));
      } else {
	if (value_size_ != src.value_size_) {
	  value_ = static_cast<void*>(realloc(value_,
					      src.value_size_ * sizeof(char)));
	}
      }
      memcpy(value_, src.value_, src.value_size_);
    } else {
      if (value_ != NULL || data_set_) {
	throw runtime_error("Data is already set value.");
      }
      if (src.value_ != NULL) {
	value_ = static_cast<void*>(calloc(src.value_size_, sizeof(char)));
	memcpy(value_, src.value_, src.value_size_);
      }
    }
    value_size_ = src.value_size_;
    value_allocated_ = true;
    data_set_ = true;
#ifdef BACKEND_KMR
    owner_ = src.owner_;
    shared_ = src.shared_;
#endif
  }

  void Data::copy_shallow(const Data& src) {
    if (value_ != NULL) {
      throw runtime_error("Data is already set value.");
    }
    value_ = src.value_;
    value_size_ = src.value_size_;
    value_allocated_ = true;
    data_set_ = true;
#ifdef BACKEND_KMR
    owner_ = src.owner_;
    shared_ = src.shared_;
#endif
  }

  void Data::copy_buf(void *val, const size_t val_siz) {
    if (value_ != NULL) {
      throw runtime_error("Data is already set value.");
    }
    value_ = static_cast<void*>(calloc(val_siz, sizeof(char)));
    memcpy(value_, val, val_siz);
    value_size_ = val_siz;
    value_allocated_ = true;
    data_set_ = true;
  }

  void Data::clear() {
    if (value_allocated_) {
      free(value_);
      value_ = NULL;
    }
    value_size_ = 0;
    value_allocated_ = false;
    data_set_ = false;
    data_file_offset_ = 0;
    data_file_size_ = 0;
#ifdef BACKEND_KMR
    owner_ = -1;
    shared_ = false;
#endif
  }

  void Data::restore_from_file_buf(char *buf) {
    if (data_file_size_ == 0) {
      throw runtime_error("Data can not be found in the file.");
    }
    copy_buf(buf + data_file_offset_, data_file_size_);
  }

  void Data::written(size_t start_pos, size_t written_siz, bool clear_data) {
    if (clear_data) {
      clear();
    }
    data_set_ = true;
    data_file_offset_ = start_pos;
    data_file_size_ = written_siz;
  }

  bool Data::updated_in_memory() {
    if (data_set_ && value_allocated_) {
      return true;
    } else {
      return false;
    }
  }

  bool Data::removed_in_memory() {
    if (!data_set_ && !value_allocated_) {
      return true;
    } else {
      return false;
    }
  }

  void Data::clear_cache() {
    if (value_allocated_) {
      free(value_);
      value_ = NULL;
    }
    value_size_ = 0;
    value_allocated_ = false;
  }

  bool Data::operator==(const Data& rhs) const {
    if (value_size_ != rhs.value_size_) {
      return false;
    }
    int cc = memcmp(value_, rhs.value_, value_size_);
    if (cc != 0) {
      return false;
    }
    return true;
  }

  DataPack::DataPack(const Key k, Data *d, bool allocate)
    : key_(k), allocated_(allocate) {
    if (allocated_) {
      data_ = new Data();
      data_->set_value(*d);
    } else {
      data_ = d;
    }
  }

  DataPack::DataPack(const DataPack &obj)
    : key_(obj.key_), allocated_(obj.allocated_) {
    if (allocated_) {
      data_ = new Data();
      data_->set_value(*(obj.data_));
    } else {
      data_ = obj.data_;
    }
  }

}
