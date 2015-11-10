#ifndef KMRNEXT_IMPL_HPP
#define KMRNEXT_IMPL_HPP

#include "../config.h"

#include "kmrnext.hpp"

#ifdef BACKEND_KMR
#include <mpi.h>
#include <kmr.h>
#endif

namespace kmrnext {

#ifdef BACKEND_KMR
  class KMRNextKmr : public KMRNext {
  public:
    KMRNextKmr();

    virtual ~KMRNextKmr();

  private:
    MPI_Comm world_comm_;
    KMR *mr_;
    int nprocs_;
    int rank_;
  };
#endif

}

#endif
