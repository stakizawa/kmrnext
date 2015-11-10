#include "../config.hpp"
#include "kmrnext.hpp"

namespace kmrnext {

  KMRNext *KMRNext::kmrnext_ = NULL;

  KMRNext* KMRNext::init(int argc, char **argv) {
    if (kmrnext_ == NULL) {
      MPI_Init(&argc, &argv);
      kmr_init();
      kmrnext_ = new KMRNext();
    }
    return kmrnext_;
  }

  void KMRNext::finalize() {
    if (kmrnext_ != NULL) {
      delete kmrnext_;
      kmr_fin();
      MPI_Finalize();
    }
  }

  DataStore* KMRNext::create_ds(size_t siz) {
    return new DataStore(siz, this);
  }

  KMRNext::KMRNext() {
    world_comm_ = MPI_COMM_WORLD;
    mr_ = kmr_create_context(world_comm_, MPI_INFO_NULL, NULL);
    MPI_Comm_size(world_comm_, &nprocs_);
    MPI_Comm_rank(world_comm_, &rank_);
  }

  KMRNext::~KMRNext() {
    kmr_free_context(mr_);
  }

}
