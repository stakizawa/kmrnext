#ifndef UTIL_HPP
#define UTIL_HPP
/// \file
/// Utility functions

#include <iostream>
#include <fstream>
#include <iomanip>
#include <ctime>
#include <sys/stat.h>
#include <unistd.h>
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

// It returns true if the file exists.
bool file_exist(string &filename) {
  struct stat st;
  int ret = stat(filename.c_str(), &st);
  if (ret == 0) {
    return true;
  } else {
    return false;
  }
}

size_t file_size(string &filename) {
  ifstream fin;
  fin.open(filename.c_str(), ios::in|ios::binary);
  fin.seekg(0, ios::end);
  size_t siz = fin.tellg();
  fin.close();
  return siz;
}

void delete_file(string &filename) {
  if (file_exist(filename)) {
    unlink(filename.c_str());
  }
}

#endif
