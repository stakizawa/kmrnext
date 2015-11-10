#include "../config.h"

#include "kmrnext.hpp"
#include "kmrnext_impl.hpp"

#ifdef BACKEND_KMR
#include <mpi.h>
#include <kmr.h>
#endif

namespace kmrnext {

  KMRNext *KMRNext::kmrnext_ = NULL;

  KMRNext* KMRNext::init(int argc, char **argv) {
    if (kmrnext_ == NULL) {
#ifdef BACKEND_KMR
      MPI_Init(&argc, &argv);
      kmr_init();
      kmrnext_ = new KMRNextKmr();
#else
      kmrnext_ = new KMRNext();
#endif
    }
    return kmrnext_;
  }

  void KMRNext::finalize() {
    if (kmrnext_ != NULL) {
      delete kmrnext_;
#ifdef BACKEND_KMR
      kmr_fin();
      MPI_Finalize();
#endif
    }
  }

  DataStore* KMRNext::create_ds(size_t siz) {
    return new DataStore(siz, this);
  }

#ifdef BACKEND_KMR
  KMRNextKmr::KMRNextKmr() {
    world_comm_ = MPI_COMM_WORLD;
    mr_ = kmr_create_context(world_comm_, MPI_INFO_NULL, NULL);
    MPI_Comm_size(world_comm_, &nprocs_);
    MPI_Comm_rank(world_comm_, &rank_);
  }

  KMRNextKmr::~KMRNextKmr() {
    kmr_free_context(mr_);
  }
#endif

}
