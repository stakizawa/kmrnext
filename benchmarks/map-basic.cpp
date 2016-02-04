/// \file
/// Benchmark program for evaluating DataStore::map() operation
///
/// It requires KMR backend.
/// It measures time for DataStore::map() by changing the following
/// parameters; # of dimensions, View type, # of data count in a DataStore
/// and size of each data.  This program accepts the following options to
/// perform the evaluation.
///    - dimension
///    - view
///    - count
///    - size
///
/// The number of processes to run this program should be two for DEBUG mode
/// and 1000 for large scale run.
#include <iostream>
#include <sstream>
#include <cstdlib>
#include <cassert>
#include <ctime>
#include "kmrnext.hpp"

using namespace std;
using namespace kmrnext;

#if DEBUG
const int kNumProcs      = 2;
const int kNumIterations = 2;
#else
const int kNumProcs      = 1000;
const int kNumIterations = 10;
#endif

void setup(int argc, int* rankp, int* nprocsp);
void eval_dimension(KMRNext* next, int rank, int nprocs);
void eval_view(KMRNext* next, int rank, int nprocs);
void eval_count(KMRNext* next, int rank, int nprocs);
void eval_size(KMRNext* next, int rank, int nprocs);
void help();


//////////////////////////////////////////////////////////////////////////////
// Main starts from here.
//////////////////////////////////////////////////////////////////////////////
int
main(int argc, char** argv)
{
  KMRNext *next = KMRNext::init(argc, argv);
  //next->enable_profile();
  int rank, nprocs;
  setup(argc, &rank, &nprocs);
  string option(argv[1]);
  if (option == "dimension") {
    eval_dimension(next, rank, nprocs);
  } else if (option == "view") {
    eval_view(next, rank, nprocs);
  } else if (option == "count") {
    eval_count(next, rank, nprocs);
  } else if (option == "size") {
    eval_size(next, rank, nprocs);
  } else {
    help();
    return 0;
  }

  KMRNext::finalize();
  return 0;
}

class Timer {
  vector<double> starts;
  vector<double> finishs;

  double gettime() {
    MPI_Barrier(MPI_COMM_WORLD);
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    return ((double) ts.tv_sec) * 10E9 + ((double) ts.tv_nsec);
  }

public:
  void start() {
    starts.push_back(gettime());
  }

  void finish() {
    finishs.push_back(gettime());
  }

  string str(bool header=true, bool header_only=false) {
    assert(starts.size() == finishs.size());
    ostringstream os;
    if (header || header_only) {
      for (size_t i = 0; i < starts.size(); i++) {
	os << i << ",";
      }
      os << "average" << endl;
      if (header_only) {
	return os.str();
      }
    }
    double time_sum = 0;
    for (size_t i = 0; i < starts.size(); i++) {
      double time = finishs.at(i) - starts.at(i);
      time_sum += time;
      os << time << ",";
    }
    os << (time_sum / (double)starts.size()) << endl;
    return os.str();
  }
};

class DataLoader : public DataStore::Loader<long> {
  size_t element_count_;

public:
  DataLoader(size_t count) : element_count_(count) {}

  int operator()(DataStore* ds, const long& num)
  {
    size_t data_count = 1;
    for (size_t i = 0; i < ds->size(); i++) {
      data_count *= ds->dim(i);
    }

    long *value = new long[element_count_];
    for (size_t i = 0; i < element_count_; i++) {
      value[i] = num;
    }
    Data data((void*)value, sizeof(long) * element_count_);
    delete[] value;

    for (size_t i = 0; i < data_count; i++) {
      Key key = ds->index_to_key(i);
      ds->add(key, data);
    }
    return 0;
  }
};

class ByteLoader : public DataStore::Loader<long> {
  size_t bytes_;

public:
  ByteLoader(size_t bytes) : bytes_(bytes) {}

  int operator()(DataStore* ds, const long& num)
  {
    size_t data_count = 1;
    for (size_t i = 0; i < ds->size(); i++) {
      data_count *= ds->dim(i);
    }

    size_t data_siz = bytes_ / sizeof(long);
    long *value = new long[data_siz];
    for (size_t i = 0; i < data_siz; i++) {
      value[i] = num;
    }
    Data data((void*)value, bytes_);
    delete[] value;

    for (size_t i = 0; i < data_count; i++) {
      Key key = ds->index_to_key(i);
      ds->add(key, data);
    }
    return 0;
  }
};

class PseudoMapper : public DataStore::Mapper {
  int world_rank_;

public:
  PseudoMapper(int rank) : world_rank_(rank) {}

  int operator()(DataStore* inds, DataStore* outds,
		 Key& key, vector<DataPack>& dps,
		 DataStore::MapEnvironment& env)
  {
#if 0
#if DEBUG
    size_t local_count = dps.size();
    size_t total_count;
    MPI_Allreduce(&local_count, &total_count, 1, MPI_LONG, MPI_SUM,
		  env.mpi_comm);
    if (world_rank_ == 0) {
      cerr << "Total data count is " << total_count << endl;
    }
#endif
#endif
    return 0;
  }
};


void
setup(int argc, int* rankp, int* nprocsp)
{
  if (argc != 2) {
    help();
    exit(1);
  }
  MPI_Comm_rank(MPI_COMM_WORLD, rankp);
  MPI_Comm_size(MPI_COMM_WORLD, nprocsp);
  if (*nprocsp < kNumProcs) {
    if (*rankp == 0) {
      cerr << "Number of processes should be bigger than " << kNumProcs
	   << endl;
    }
    KMRNext::finalize();
    exit(1);
  }
  if (*nprocsp != kNumProcs) {
    if (*rankp == 0) {
      cerr << "Number of runnning processes is bigger than the required "
	   << "number of processes." << endl;
      cerr << (*nprocsp - kNumProcs) << " processes are not used." << endl;
    }
  }
}

static void
eval_dimension_once(string description, KMRNext* next, int rank, int ndim,
		    size_t* dim_ary, bool* view_ary,
		    vector<long>& datalist, DataLoader& loader)
{
  DataStore *ds = next->create_ds(ndim);
  ds->set(dim_ary);
  ds->load_array(datalist, loader);
  View view(ndim);
  view.set(view_ary);
  Timer timer;
  for (size_t i = 0; i < kNumIterations; i++) {
    PseudoMapper mapper(rank);
    timer.start();
    ds->map(NULL, mapper, view);
    timer.finish();
  }
  delete ds;
  if (rank == 0) {
    cout << description << endl;
    cout << timer.str() << endl;
  }
}

void
eval_dimension(KMRNext* next, int rank, int nprocs)
{
  if (rank == 0) {
    cout << "Dimension Test" << endl;
  }

#if DEBUG
  size_t ary2[2] = {kNumProcs, 10000};
  size_t ary3[3] = {kNumProcs, 10, 1000};
  size_t ary4[4] = {kNumProcs, 10, 10, 100};
  size_t ary5[5] = {kNumProcs, 10, 10, 10, 10};
  size_t ary6[6] = {kNumProcs, 10, 10, 10, 5, 2};
  size_t ary7[7] = {kNumProcs, 10, 10, 5, 2, 5, 2};
  size_t ary8[8] = {kNumProcs, 10, 5, 2, 5, 2, 5, 2};
#else
  size_t ary2[2] = {kNumProcs, 10000};
  size_t ary3[3] = {kNumProcs, 10, 1000};
  size_t ary4[4] = {kNumProcs, 10, 10, 100};
  size_t ary5[5] = {kNumProcs, 10, 10, 10, 10};
  size_t ary6[6] = {kNumProcs, 10, 10, 10, 5, 2};
  size_t ary7[7] = {kNumProcs, 10, 10, 5, 2, 5, 2};
  size_t ary8[8] = {kNumProcs, 10, 5, 2, 5, 2, 5, 2};
#endif

  bool   aryv2[2] = {true, false};
  bool   aryv3[3] = {true, false, false};
  bool   aryv4[4] = {true, false, false, false};
  bool   aryv5[5] = {true, false, false, false, false};
  bool   aryv6[6] = {true, false, false, false, false, false};
  bool   aryv7[7] = {true, false, false, false, false, false, false};
  bool   aryv8[8] = {true, false, false, false, false, false, false, false};

  vector<long> datalist;
  for (size_t i = 0; i < kNumProcs; i++) {
    datalist.push_back((long)(i+1));
  }
  DataLoader loader(1);

  eval_dimension_once(string("2D"),
		      next, rank, 2, ary2, aryv2, datalist, loader);
  eval_dimension_once(string("3D"),
		      next, rank, 3, ary3, aryv3, datalist, loader);
  eval_dimension_once(string("4D"),
		      next, rank, 4, ary4, aryv4, datalist, loader);
  eval_dimension_once(string("5D"),
		      next, rank, 5, ary5, aryv5, datalist, loader);
  eval_dimension_once(string("6D"),
		      next, rank, 6, ary6, aryv6, datalist, loader);
  eval_dimension_once(string("7D"),
		      next, rank, 7, ary7, aryv7, datalist, loader);
  eval_dimension_once(string("8D"),
		      next, rank, 8, ary8, aryv8, datalist, loader);
}

static void
eval_view_once(string description, int rank, DataStore* ds,
	       bool* view_ary)
{
  View view(ds->size());
  view.set(view_ary);
  Timer timer;
  for (size_t i = 0; i < kNumIterations; i++) {
    PseudoMapper mapper(rank);
    timer.start();
    ds->map(NULL, mapper, view);
    timer.finish();
  }
  if (rank == 0) {
    cout << description << endl;
    cout << timer.str() << endl;
  }
}

void
eval_view(KMRNext* next, int rank, int nprocs)
{
  if (rank == 0) {
    cout << "View Test" << endl;
  }

#if DEBUG
  size_t ary3[3] = {kNumProcs, 100, 100};
#else
  size_t ary3[3] = {kNumProcs, 100, 100};
#endif

  bool fff[3] = {false, false, false};
  bool fft[3] = {false, false, true };
  bool ftf[3] = {false, true,  false};
  bool ftt[3] = {false, true,  true };
  bool tff[3] = {true,  false, false};
  bool tft[3] = {true,  false, true };
  bool ttf[3] = {true,  true,  false};
  bool ttt[3] = {true,  true,  true };

  vector<long> datalist;
  for (size_t i = 0; i < kNumProcs; i++) {
    datalist.push_back((long)(i+1));
  }
  DataLoader loader(1);
  DataStore *ds = next->create_ds(3);
  ds->set(ary3);
  ds->load_array(datalist, loader);

  eval_view_once(string("FFF"), rank, ds, fff);
  eval_view_once(string("FFT"), rank, ds, fft);
  eval_view_once(string("FTF"), rank, ds, ftf);
  eval_view_once(string("FTT"), rank, ds, ftt);
  eval_view_once(string("TFF"), rank, ds, tff);
  eval_view_once(string("TFT"), rank, ds, tft);
  eval_view_once(string("TTF"), rank, ds, ttf);
  eval_view_once(string("TTT"), rank, ds, ttt);

  delete ds;
}

static void
eval_count_once(string description, KMRNext* next, int rank, size_t* dim_ary,
		vector<long>& datalist, DataLoader& loader)
{
  DataStore *ds = next->create_ds(2);
  ds->set(dim_ary);
  ds->load_array(datalist, loader);

  bool ary_ff[2] = {false, false};
  View vff(2);
  vff.set(ary_ff);
  bool ary_ft[2] = {false, true};
  View vft(2);
  vft.set(ary_ft);
  bool ary_tf[2] = {true,  false};
  View vtf(2);
  vtf.set(ary_tf);
  bool ary_tt[2] = {true,  true};
  View vtt(2);
  vtt.set(ary_tt);

  Timer timer_ff;
  for (size_t i = 0; i < kNumIterations; i++) {
    PseudoMapper mapper(rank);
    timer_ff.start();
    ds->map(NULL, mapper, vff);
    timer_ff.finish();
  }
  Timer timer_ft;
  for (size_t i = 0; i < kNumIterations; i++) {
    PseudoMapper mapper(rank);
    timer_ft.start();
    ds->map(NULL, mapper, vft);
    timer_ft.finish();
  }
  Timer timer_tf;
  for (size_t i = 0; i < kNumIterations; i++) {
    PseudoMapper mapper(rank);
    timer_tf.start();
    ds->map(NULL, mapper, vtf);
    timer_tf.finish();
  }
  Timer timer_tt;
  for (size_t i = 0; i < kNumIterations; i++) {
    PseudoMapper mapper(rank);
    timer_tt.start();
    ds->map(NULL, mapper, vtt);
    timer_tt.finish();
  }
  delete ds;
  if (rank == 0) {
    cout << description << endl;
    cout << "  ," << timer_ff.str(true, true);
    cout << "FF," << timer_ff.str(false);
    cout << "FT," << timer_tf.str(false);
    cout << "TF," << timer_ft.str(false);
    cout << "TT," << timer_tt.str(false) << endl;
  }
}

void
eval_count(KMRNext* next, int rank, int nprocs)
{
  if (rank == 0) {
    cout << "Data Count Test" << endl;
  }

#if DEBUG
  size_t ary1[2] = {kNumProcs, 10};
  size_t ary2[2] = {kNumProcs, 20};
  size_t ary3[2] = {kNumProcs, 40};
  size_t ary4[2] = {kNumProcs, 80};
  size_t ary5[2] = {kNumProcs, 160};
  size_t ary6[2] = {kNumProcs, 320};
  size_t ary7[2] = {kNumProcs, 640};
  size_t ary8[2] = {kNumProcs, 1280};
#else
  size_t ary1[2] = {kNumProcs, 10};
  size_t ary2[2] = {kNumProcs, 20};
  size_t ary3[2] = {kNumProcs, 40};
  size_t ary4[2] = {kNumProcs, 80};
  size_t ary5[2] = {kNumProcs, 160};
  size_t ary6[2] = {kNumProcs, 320};
  size_t ary7[2] = {kNumProcs, 640};
  size_t ary8[2] = {kNumProcs, 1280};
#endif

  vector<long> datalist;
  for (size_t i = 0; i < kNumProcs; i++) {
    datalist.push_back((long)(i+1));
  }
  DataLoader loader(1);

  eval_count_once(string("10"), next, rank, ary1, datalist, loader);
  eval_count_once(string("20"), next, rank, ary2, datalist, loader);
  eval_count_once(string("30"), next, rank, ary3, datalist, loader);
  eval_count_once(string("40"), next, rank, ary4, datalist, loader);
  eval_count_once(string("50"), next, rank, ary5, datalist, loader);
  eval_count_once(string("60"), next, rank, ary6, datalist, loader);
  eval_count_once(string("70"), next, rank, ary7, datalist, loader);
  eval_count_once(string("80"), next, rank, ary8, datalist, loader);
}

static void
eval_size_once(string description, KMRNext* next, int rank, size_t data_size,
	       vector<long>& datalist)
{
  DataStore *ds = next->create_ds(2);
  size_t dim_ary[2] = {kNumProcs, 100};
  ds->set(dim_ary);
  ByteLoader loader(data_size);
  ds->load_array(datalist, loader);

  bool ary_ff[2] = {false, false};
  View vff(2);
  vff.set(ary_ff);
  bool ary_ft[2] = {false, true};
  View vft(2);
  vft.set(ary_ft);
  bool ary_tf[2] = {true,  false};
  View vtf(2);
  vtf.set(ary_tf);
  bool ary_tt[2] = {true,  true};
  View vtt(2);
  vtt.set(ary_tt);

  Timer timer_ff;
  for (size_t i = 0; i < kNumIterations; i++) {
    PseudoMapper mapper(rank);
    timer_ff.start();
    ds->map(NULL, mapper, vff);
    timer_ff.finish();
  }
  Timer timer_ft;
  for (size_t i = 0; i < kNumIterations; i++) {
    PseudoMapper mapper(rank);
    timer_ft.start();
    ds->map(NULL, mapper, vft);
    timer_ft.finish();
  }
  Timer timer_tf;
  for (size_t i = 0; i < kNumIterations; i++) {
    PseudoMapper mapper(rank);
    timer_tf.start();
    ds->map(NULL, mapper, vtf);
    timer_tf.finish();
  }
  Timer timer_tt;
  for (size_t i = 0; i < kNumIterations; i++) {
    PseudoMapper mapper(rank);
    timer_tt.start();
    ds->map(NULL, mapper, vtt);
    timer_tt.finish();
  }
  delete ds;
  if (rank == 0) {
    cout << description << endl;
    cout << "  ," << timer_ff.str(true, true);
    cout << "FF," << timer_ff.str(false);
    cout << "FT," << timer_tf.str(false);
    cout << "TF," << timer_ft.str(false);
    cout << "TT," << timer_tt.str(false) << endl;
  }
}

void
eval_size(KMRNext* next, int rank, int nprocs)
{
  if (rank == 0) {
    cout << "Data Size Test" << endl;
  }

  size_t sizes[11] = {1024, 2048, 4096, 8192, 16384, 32768, 65536,
		      131072, 262144, 524288, 1048576};

  vector<long> datalist;
  for (size_t i = 0; i < kNumProcs; i++) {
    datalist.push_back((long)(i+1));
  }

  eval_size_once(string("1K"), next, rank, sizes[0], datalist);
  eval_size_once(string("2K"), next, rank, sizes[1], datalist);
  eval_size_once(string("4K"), next, rank, sizes[2], datalist);
  eval_size_once(string("8K"), next, rank, sizes[3], datalist);
#if not DEBUG
  eval_size_once(string("16K"), next, rank, sizes[4], datalist);
  eval_size_once(string("32K"), next, rank, sizes[5], datalist);
  eval_size_once(string("64K"), next, rank, sizes[6], datalist);
  eval_size_once(string("128K"), next, rank, sizes[7], datalist);
  eval_size_once(string("256K"), next, rank, sizes[8], datalist);
  eval_size_once(string("512K"), next, rank, sizes[9], datalist);
  eval_size_once(string("1024K"), next, rank, sizes[10], datalist);
#endif
}

void
help()
{
  cout << "Options: " << endl;
  cout << "  dimension" << endl;
  cout << "      Change number of dimensions of a DataStore" << endl;
  cout << "  view" << endl;
  cout << "      Change view type" << endl;
  cout << "  count" << endl;
  cout << "      Change count of data in a DataStore" << endl;
  cout << "  size" << endl;
  cout << "      Change size of data in a DataStore" << endl;
  cout << endl;
  KMRNext::finalize();
}
