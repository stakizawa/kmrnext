#include "../config.hpp"
#include "kmrnext.hpp"

namespace kmrnext {

  KMRNext *KMRNext::kmrnext_ = NULL;

  DataStore* KMRNext::create_ds(size_t siz) {
    return new DataStore(siz, this);
  }

}

#ifdef BACKEND_SERIAL
#include "kmrnext_serial.cpp"
#elif defined BACKEND_KMR
#include "kmrnext_kmr.cpp"
#endif
