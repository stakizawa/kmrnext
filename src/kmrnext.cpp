#include "../config.hpp"
#include "kmrnext.hpp"

#ifdef BACKEND_SERIAL
#include "kmrnext_serial.cpp"
#elif defined BACKEND_KMR
#include "kmrnext_kmr.cpp"
#endif
