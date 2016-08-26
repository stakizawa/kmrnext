#include "../config.hpp"
#include "kmrnext.hpp"

namespace kmrnext {

  KMRNext* KMRNext::init(int argc, char **argv) {
    return KMRNext::init();
  }

  KMRNext* KMRNext::init() {
    if (kmrnext_ == NULL) {
      kmrnext_ = new KMRNext();
    }
    DataStore::initialize(kmrnext_);
    return kmrnext_;
  }

  void KMRNext::finalize() {
    DataStore::finalize();
    if (kmrnext_ != NULL) {
      delete kmrnext_;
      kmrnext_ = NULL;
    }
  }

  KMRNext::KMRNext() : profile_(false), iomode_(KMRNext::Memory) {}

  KMRNext::~KMRNext() {}

  void KMRNext::enable_profile() {
    profile_ = true;
  }

  void KMRNext::disable_profile() {
    profile_ = false;
  }

}
