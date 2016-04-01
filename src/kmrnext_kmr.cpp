#ifdef _OPENMP
#include <omp.h>
#endif
#include "../config.hpp"
#include "kmrnext.hpp"

namespace kmrnext {

  KMRNext* KMRNext::init(int argc, char **argv) {
    if (kmrnext_ == NULL) {
#ifdef _OPENMP
      int thlv;
      MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &thlv);
      if (thlv < MPI_THREAD_SERIALIZED) {
	cerr << "Warning:  This MPI implementation provides insufficient"
	     << " threading support." << endl;
	omp_set_num_threads(1);
      }
#else
      MPI_Init(&argc, &argv);
#endif
      kmr_init();
      kmrnext_ = new KMRNext();
    }
    DataStore::initialize();
    return kmrnext_;
  }

  void KMRNext::finalize() {
    DataStore::finalize();
    if (kmrnext_ != NULL) {
      delete kmrnext_;
      kmr_fin();
      MPI_Finalize();
    }
  }

  KMRNext::KMRNext() {
    world_comm_ = MPI_COMM_WORLD;
    mr_ = kmr_create_context(world_comm_, MPI_INFO_NULL, NULL);
    MPI_Comm_size(world_comm_, &nprocs_);
    MPI_Comm_rank(world_comm_, &rank_);
    profile_ = false;
  }

  KMRNext::~KMRNext() {
    kmr_free_context(mr_);
  }

  void KMRNext::enable_profile() {
    profile_ = true;
    mr_->verbosity = 5; // larger than 7 outputs more verbosely
    mr_->trace_map_mp = 1;
  }

  void KMRNext::disable_profile() {
    profile_ = false;
    mr_->verbosity = 5; // larger than 7 outputs more verbosely
    mr_->trace_map_mp = 0;
  }

}
