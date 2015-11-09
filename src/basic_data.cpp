#include <cstdlib>
#include <cstring>
#include "kmrnext.hpp"

namespace kmrnext {

  Data::~Data() {
    if (value_allocated_) {
      free(value_);
    }
  }

  void Data::copy_deep(const Data& src) {
    if (value_ != NULL) {
      throw runtime_error("Data is already set value.");
    }
    value_ = static_cast<void*>(calloc(src.value_size_, sizeof(void*)));
    memcpy(value_, src.value_, src.value_size_);
    value_size_ = src.value_size_;
    value_allocated_ = true;
  }

}
