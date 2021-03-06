#include "../config.hpp"
#include "kmrnext.hpp"
#ifdef BACKEND_SERIAL
#include "data_store_serial.hpp"
#elif defined BACKEND_KMR
#include "data_store_kmr.hpp"
#endif

namespace kmrnext {

  KMRNext* KMRNext::kmrnext_ = NULL;

  DataStore* KMRNext::create_ds(size_t siz) {
    if (iomode_ == KMRNext::Memory) {
      return new InMemoryDataStore(siz, this);
    } else if (iomode_ == KMRNext::File) {
      return new SimpleFileDataStore(siz, this);
    } else {
      throw runtime_error("Unknown media type.");
    }
  }

}

#ifdef BACKEND_SERIAL
#include "kmrnext_serial.cpp"
#elif defined BACKEND_KMR
#include "kmrnext_kmr.cpp"
#endif
