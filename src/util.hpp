#ifndef UTIL_HPP
#define UTIL_HPP
/// \file
/// Utility functions

#include <iostream>
#include <iomanip>
#include <ctime>
#include "kmrnext.hpp"

using namespace std;

#ifdef BACKEND_KMR
// It returns current time in nano seconds after synchronizing between
// processes in the specified MPI communicator.
double gettime(MPI_Comm comm) {
  MPI_Barrier(comm);
  struct timespec ts;
  clock_gettime(CLOCK_REALTIME, &ts);
  return ((double)ts.tv_sec) * 10E9 + ((double)ts.tv_nsec);
}
#endif

#ifdef BACKEND_SERIAL
void profile_out(string message) {
  cerr << ";;KMRNEXT: " << message << endl;
}
#elif defined BACKEND_KMR
void profile_out(kmrnext::KMRNext *kmrnext_, string message) {
  cerr << ";;KMRNEXT [" << setfill('0') << setw(5) << kmrnext_->rank() << "] "
       << message << endl;
}
#endif

#endif
