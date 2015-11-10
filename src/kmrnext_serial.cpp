#include "../config.hpp"
#include "kmrnext.hpp"

namespace kmrnext {

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

  KMRNext::KMRNext() {}

  KMRNext::~KMRNext() {}

}
