#include "../config.hpp"

#include <cstdlib>
#include <cstring>
#include "kmrnext.hpp"

namespace kmrnext {

  string View::to_string() const {
    ostringstream os;
    os << '<';
    for (size_t i = 0; i < size_; i++) {
      if (value_[i] == SplitAll) {
	os << "All";
      } else if (value_[i] == SplitNone) {
	os << "None";
      } else {
	os << value_[i];
      }
      if (i < size_ - 1) {
	os << ',';
      }
    }
    os << '>';
    return os.str();
  }

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
    void* tmp = value_;
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

}
