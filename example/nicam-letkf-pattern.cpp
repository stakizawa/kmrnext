#include <iostream>
#include <sstream>
#include <cassert>
#include "kmrnext.hpp"

using namespace std;
using namespace kmrnext;

int rank = 0;

const size_t kNumIteration = 10;

const size_t kDimEnsembleData = 3;
//const size_t kDimRegionData   = 2;
const size_t kDimLatticeData  = 1;

const size_t kNumEnsemble = 2;
const size_t kNumRegion   = 10;
#if DEBUG
const size_t kNumLattice  = 10;
#else
const size_t kNumLattice  = 1156;
#endif

const size_t kEnsembleDataDimSizes[kDimEnsembleData] =
  {kNumEnsemble, kNumRegion, kNumLattice};

const bool kPrint = true;

void load_data(DataStore* ds);
void run_nicam(DataStore* inds, DataStore* outds);
void run_letkf(DataStore* inds, DataStore* outds);
void print_line(string& str);
void print_line(ostringstream& os);

class DataPrinter : public DataPack::Dumper {
public:
  string operator()(DataPack& dp)
  {
    ostringstream os;
    os << dp.key().to_string() << " : " << *(int*)(dp.data()->value()) << endl;
    return os.str();
  }
};

//////////////////////////////////////////////////////////////////////////////
// Main starts from here.
//////////////////////////////////////////////////////////////////////////////
int
main(int argc, char **argv)
{
  KMRNext *next = KMRNext::init(argc, argv);
#ifdef BACKEND_KMR
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  DataStore* ds0 = next->create_ds(kDimEnsembleData);
  ds0->set(kEnsembleDataDimSizes);
  load_data(ds0);
  DataPrinter dp;
#if DEBUG
  string str0 = ds0->dump(dp);
  print_line(str0);
#endif

  for (size_t i = 0; i < kNumIteration; i++) {
    ostringstream os0;
    os0 << "Iteration[" << i << "] starts.";
    print_line(os0);

    // run pseudo-NICAM
    DataStore* ds1 = next->create_ds(kDimEnsembleData);
    ds1->set(kEnsembleDataDimSizes);
    run_nicam(ds0, ds1);
    delete ds0;
#if DEBUG
    string str1 = ds1->dump(dp);
    print_line(str1);
#endif

    // run pseudo-LETKF
    ds0 = next->create_ds(kDimEnsembleData);
    ds0->set(kEnsembleDataDimSizes);
    run_letkf(ds1, ds0);
    delete ds1;
#if DEBUG
    string str2 = ds0->dump(dp);
    print_line(str2);
#endif

    ostringstream os1;
    os1 << "Iteration[" << i << "] ends.";
    print_line(os1);
  }

  KMRNext::finalize();
  return 0;
}


class DataLoader : public DataStore::Loader<int> {
public:
  int operator()(DataStore* ds, const int& num)
  {
    Key key(kDimLatticeData);
    Data data((void*)&num, sizeof(int));
    for (size_t i = 0; i < kNumLattice; i++) {
      key.set_dim(0, i);
      ds->add(key, data);
    }
    return 0;
  }
};

void load_data(DataStore* ds)
{
  vector<int> data_srcs;
  for (size_t i = 0; i < kNumEnsemble; i++) {
    for (size_t j = 0; j < kNumRegion; j++) {
      data_srcs.push_back((int)(j+1));
    }
  }
  DataLoader dl;
  ds->load_array(data_srcs, dl);
}

class PseudoNICAM : public DataStore::Mapper {
public:
  int operator()(DataStore* inds, DataStore* outds,
		 Key& key, vector<DataPack>& dps,
		 DataStore::MapEnvironment& env)
  {
#ifdef BACKEND_SERIAL
    assert(dps.size() == (size_t)(kNumRegion * kNumLattice));
#elif defined BACKEND_KMR
    size_t local_count = dps.size();
    size_t total_count;
    MPI_Allreduce(&local_count, &total_count, 1, MPI_LONG, MPI_SUM,
		  env.mpi_comm);
    assert(total_count == (size_t)(kNumRegion * kNumLattice));
#endif

    for (vector<DataPack>::iterator itr = dps.begin(); itr != dps.end();
	 itr++) {
      int newval = *(int*)itr->data()->value() + 1;
      Data data(&newval, itr->data()->size());
      outds->add(itr->key(), data);
    }
    return 0;
  }
};

void run_nicam(DataStore* inds, DataStore* outds)
{
  PseudoNICAM mapper;
  View view(kDimEnsembleData);
  bool view_flag[3] = {true, false, false};
  view.set(view_flag);
  inds->map(outds, mapper, view);
}

class PseudoLETKF : public DataStore::Mapper {
public:
  int operator()(DataStore* inds, DataStore* outds,
		 Key& key, vector<DataPack>& dps,
		 DataStore::MapEnvironment& env)
  {
#ifdef BACKEND_SERIAL
    assert(dps.size() == (size_t)kNumEnsemble);
#elif defined BACKEND_KMR
    size_t local_count = dps.size();
    size_t total_count;
    MPI_Allreduce(&local_count, &total_count, 1, MPI_LONG, MPI_SUM,
		  env.mpi_comm);
    assert(total_count == (size_t)kNumEnsemble);
#endif

    for (vector<DataPack>::iterator itr = dps.begin(); itr != dps.end();
	 itr++) {
      int newval = *(int*)itr->data()->value() - 1;
      Data data(&newval, itr->data()->size());
      outds->add(itr->key(), data);
    }
    return 0;
  }
};

void run_letkf(DataStore* inds, DataStore* outds)
{
  PseudoLETKF mapper;
  View view(kDimEnsembleData);
  bool view_flag[3] = {false, true, true};
  view.set(view_flag);
  inds->map(outds, mapper, view);
}

void print_line(string& str) {
  if (kPrint && (rank == 0)) {
    cout << str << endl;
  }
}

void print_line(ostringstream& os) {
  string s = os.str();
  print_line(s);
}
