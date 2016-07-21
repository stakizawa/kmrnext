#include "../config.hpp"

#include <cstdlib>
#include <cstring>
#include "kmrnext.hpp"

namespace kmrnext {

  Data::Data(const Data& obj)
    : value_(obj.value_), value_size_(obj.value_size_),
      value_allocated_(false) {}

  Data::Data(const Data* obj) : value_allocated_(false) {
    if (obj != NULL) {
      value_ = obj->value_;
      value_size_ = obj->value_size_;
    } else {
      value_ = NULL;
      value_size_ = 0;
    }
  }

  Data::~Data() {
    if (value_allocated_) {
      free(value_);
    }
  }

  void Data::allocate() {
    if (value_size_ == 0) {
      return;
    }
    if (value_allocated_) {
      throw runtime_error("Already allocated.");
    }
    void *tmp = value_;
    value_ = static_cast<void*>(calloc(value_size_, sizeof(char)));
    memcpy(value_, tmp, value_size_);
    value_allocated_ = true;
  }

  bool Data::operator==(const Data& rhs) const {
    if (value_size_ != rhs.value_size_) {
      return false;
    }
    if (value_ != rhs.value_) {
      int cc = memcmp(value_, rhs.value_, value_size_);
      if (cc != 0) {
	return false;
      }
    }
    return true;
  }

  DataPack::DataPack(const Key& k, const Data* d, bool allocate)
    : key_(k), data_(d), allocated_(allocate) {
    if (allocated_) {
      data_.allocate();
    }
  }

  DataPack::DataPack(const Key& k, const Data& d, bool allocate)
    : key_(k), data_(d), allocated_(allocate) {
    if (allocated_) {
      data_.allocate();
    }
  }

  DataPack::DataPack(const DataPack &obj)
    : key_(obj.key_), data_(obj.data_), allocated_(obj.allocated_) {
    if (allocated_) {
      data_.allocate();
    }
  }

  DataElement::DataElement()
    : data_(NULL), data_set_(false),
      data_updated_(false), data_file_offset_(0), data_file_size_(0)
#ifdef BACKEND_KMR
      owner_(-1), shared_(false)
#endif
  {}

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

#if 0   // TODO
#ifdef BACKEND_KMR
    owner_ = dat.owner_;
    shared_ = dat.shared_;
#endif
#endif
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
