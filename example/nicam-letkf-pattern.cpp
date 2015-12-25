#include <iostream>
#include <sstream>
#include <cassert>
#include "kmrnext.hpp"

using namespace std;
using namespace kmrnext;

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
bool to_print();

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

  DataStore* ds0 = next->create_ds(kDimEnsembleData);
  ds0->set(kEnsembleDataDimSizes);
  load_data(ds0);
  DataPrinter dp;
#if DEBUG
  string str0 = ds0->dump(dp);
  if (to_print()) {
    cout << str0 << endl;
  }
#endif

  for (size_t i = 0; i < kNumIteration; i++) {
    if (to_print()) {
      cout << "Iteration[" << i << "] starts." << endl;
    }

    // run pseudo-NICAM
    DataStore* ds1 = next->create_ds(kDimEnsembleData);
    ds1->set(kEnsembleDataDimSizes);
    run_nicam(ds0, ds1);
    delete ds0;
#if DEBUG
    string str1 = ds1->dump(dp);
    if (to_print()) {
      cout << str1 << endl;
    }
#endif

    // run pseudo-LETKF
    ds0 = next->create_ds(kDimEnsembleData);
    ds0->set(kEnsembleDataDimSizes);
    run_letkf(ds1, ds0);
    delete ds1;
#if DEBUG
    string str2 = ds0->dump(dp);
    if (to_print()) {
      cout << str2 << endl;
    }
#endif

    if (to_print()) {
      cout << "Iteration[" << i << "] ends." << endl;
    }
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
    assert(dps.size() == (size_t)(kNumRegion * kNumLattice));

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
    assert(dps.size() == (size_t)kNumEnsemble);

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

bool to_print()
{
#ifdef BACKEND_SERIAL
  return kPrint;
#elif defined BACKEND_KMR
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  return kPrint && (rank == 0);
#endif
}
