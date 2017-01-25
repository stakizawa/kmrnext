#ifndef UTIL_HPP
#define UTIL_HPP
/// \file
/// Utility functions

#include "kmrnext.hpp"

namespace kmrnext {
#ifdef BACKEND_KMR
  // It returns current time in nano seconds after synchronizing between
  // processes in the specified MPI communicator.
  double gettime(MPI_Comm comm);
#endif

  // It write message to the standard error.
#ifdef BACKEND_SERIAL
  void profile_out(std::string message);
#elif defined BACKEND_KMR
  void profile_out(kmrnext::KMRNext* kmrnext_, std::string message);
#endif

  // It returns true if the file exists.
  bool file_exist(std::string &filename);

  // It returns size of the file.
  size_t file_size(std::string &filename);

  // It deletes the file.
  void delete_file(std::string &filename);

  // It serializes a string.
  void serialize(const string& str, char** buf, size_t* buf_siz);

  // It deserializes a string.
  void deserialize(char* buf, size_t buf_siz, string** str);

  // It serializes an integer.
  void serialize(const long& val, char** buf, size_t* buf_siz);

  // It deserializes an integer.
  void deserialize(char* buf, size_t buf_siz, long** val);
}

#endif
