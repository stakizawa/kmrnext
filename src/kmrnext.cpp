#include "kmrnext.hpp"

namespace kmrnext {

  KMRNext *KMRNext::kmrnext_ = NULL;

  KMRNext* KMRNext::init(int argc, char **argv) {
    if (kmrnext_ == NULL) {
      kmrnext_ = new KMRNext();
    }
    return kmrnext_;
  }

  void KMRNext::finalize() {
    if (kmrnext_ != NULL) {
      delete kmrnext_;
    }
  }

  DataStore* KMRNext::create_ds(size_t siz) {
    return new DataStore(siz);
  }

}
