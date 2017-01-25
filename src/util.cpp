#include <iostream>
#include <fstream>
#include <sys/stat.h>
#include <unistd.h>
#ifdef BACKEND_KMR
#include <iomanip>
#include <ctime>
#endif

#include "kmrnext.hpp"
#include "util.hpp"

namespace kmrnext {
  using namespace std;

#ifdef BACKEND_KMR
  double gettime(MPI_Comm comm) {
    MPI_Barrier(comm);
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return (static_cast<double>(ts.tv_sec) * 1E9 +
	    static_cast<double>(ts.tv_nsec));
  }
#endif

#ifdef BACKEND_SERIAL
  void profile_out(string message) {
    cerr << ";;KMRNEXT: " << message << endl;
  }
#elif defined BACKEND_KMR
  void profile_out(kmrnext::KMRNext* kmrnext_, string message) {
    cerr << ";;KMRNEXT [" << setfill('0') << setw(5) << kmrnext_->rank()
	 << "] " << message << endl;
  }
#endif

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
    streampos siz = fin.tellg();
    fin.close();
    return static_cast<size_t>(siz);
  }

  void delete_file(string &filename) {
    if (file_exist(filename)) {
      unlink(filename.c_str());
    }
  }

  void serialize(const string& str, char** buf, size_t* buf_siz) {
    *buf = const_cast<char*>(str.c_str());
    *buf_siz = str.size() + 1; // +1 for '\0'
  }

  void deserialize(char* buf, size_t buf_siz, string** str) {
    *str = new string(buf);
  }

  void serialize(const long& val, char** buf, size_t* buf_siz) {
    long* v = const_cast<long*>(&val);
    *buf = reinterpret_cast<char*>(v);
    *buf_siz = sizeof(long);
  }

  void deserialize(char* buf, size_t buf_siz, long** val) {
    *val = new long;
    **val = *reinterpret_cast<long*>(buf);
  }

}
