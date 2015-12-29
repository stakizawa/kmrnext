#include "../config.hpp"

#include <cstdlib>
#include <cstring>
#include "kmrnext.hpp"

namespace kmrnext {

#ifdef BACKEND_SERIAL
  Data::Data() :
    value_(NULL), value_size_(0), value_allocated_(false) {}

  Data::Data(void *val, const size_t val_siz) :
    value_(val), value_size_(val_siz), value_allocated_(false) {}
#elif defined BACKEND_KMR
  Data::Data() :
    value_(NULL), value_size_(0), value_allocated_(false),
    owner_(-1), shared_(false) {}

  Data::Data(void *val, const size_t val_siz) :
    value_(val), value_size_(val_siz), value_allocated_(false),
    owner_(-1), shared_(false) {}
#endif

  Data::~Data() {
    if (value_allocated_) {
      free(value_);
    }
  }

  void Data::copy_deep(const Data& src) {
    if (value_ != NULL) {
      throw runtime_error("Data is already set value.");
    }
    value_ = static_cast<void*>(calloc(src.value_size_, sizeof(char)));
    memcpy(value_, src.value_, src.value_size_);
    value_size_ = src.value_size_;
    value_allocated_ = true;
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
#ifdef BACKEND_KMR
    owner_ = src.owner_;
    shared_ = src.shared_;
#endif
  }

}
