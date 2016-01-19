#ifndef UTIL_HPP
#define UTIL_HPP
/// \file
/// Utility functions

#include <iostream>
#include <ctime>
#include "kmrnext.hpp"

using namespace std;

#if 0
/// It returns current time in nano seconds.
double gettime() {
  struct timespec ts;
  clock_gettime(CLOCK_REALTIME, &ts);
  cout << (double)ts.tv_sec << endl;
  return ((double)ts.tv_sec) * 10E9 + ((double)ts.tv_nsec);
}
#endif

void profile_out(string message) {
  cerr << ";;KMRNEXT: " << message << endl;
}

#endif
