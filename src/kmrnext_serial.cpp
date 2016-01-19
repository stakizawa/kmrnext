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

  KMRNext::KMRNext() : profile_(false) {}

  KMRNext::~KMRNext() {}

  void KMRNext::enable_profile() {
    profile_ = true;
  }

  void KMRNext::disable_profile() {
    profile_ = false;
  }

}
