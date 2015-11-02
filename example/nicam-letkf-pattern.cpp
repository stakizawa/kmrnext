#include <iostream>
#include <sstream>
#include <cassert>
#include "kmrnext.hpp"

using namespace std;
using namespace Next;

const int NumIteration = 10;

const int DimEnsembleData = 3;
const int DimRegionData   = 2;
const int DimLatticeData  = 1;

const int NumEnsemble = 2;
const int NumRegion   = 10;
#if DEBUG
const int NumLattice  = 10;
#else
const int NumLattice  = 1156;
#endif

const size_t EnsembleDataDimSizes[DimEnsembleData] =
  {NumEnsemble, NumRegion, NumLattice};

void load_data(DataStore* ds);
void run_nicam(DataStore* inds, DataStore* outds);
void run_letkf(DataStore* inds, DataStore* outds);

class DataPrinter : public DataPack::Dumper {
public:
  string operator()(DataPack& dp)
  {
    ostringstream os;
    os << dp.key.to_string() << " : " << *(int*)(dp.data->value()) << endl;
    return os.str();
  }
};

//////////////////////////////////////////////////////////////////////////////
// Main starts from here.
//////////////////////////////////////////////////////////////////////////////
int
main(int argc, char **argv)
{
  DataStore* ds0 = new DataStore(DimEnsembleData);
  ds0->set(EnsembleDataDimSizes);
  load_data(ds0);
  DataPrinter dp;
#if DEBUG
  cout << ds0->dump(dp) << endl;
#endif

  for (int i = 0; i < NumIteration; i++) {
    // run pseudo-NICAM
    DataStore* ds1 = new DataStore(DimEnsembleData);
    ds1->set(EnsembleDataDimSizes);
    run_nicam(ds0, ds1);
    delete ds0;
#if DEBUG
    cout << ds1->dump(dp) << endl;
#endif

    // run pseudo-LETKF
    ds0 = new DataStore(DimEnsembleData);
    ds0->set(EnsembleDataDimSizes);
    run_letkf(ds1, ds0);
    delete ds1;
#if DEBUG
    cout << ds0->dump(dp) << endl;
#endif
  }

  return 0;
}


class DataLoader : public DataStore::Loader<int> {
public:
  int operator()(DataStore* ds, const int& num)
  {
    Key key(DimLatticeData);
    Data data((void*)&num, sizeof(int));
    for (int i = 0; i < NumLattice; i++) {
      key.set_dim(0, i);
      ds->add(key, data);
    }
    return 0;
  }
};

void load_data(DataStore* ds)
{
  vector<int> data_srcs;
  for (int i = 0; i < NumEnsemble; i++) {
    for (int j = 0; j < NumRegion; j++) {
      data_srcs.push_back(j+1);
    }
  }
  DataLoader dl;
  ds->load_array(data_srcs, dl);
}

class PseudoNICAM : public DataStore::Mapper {
public:
  int operator()(DataStore* inds, DataStore* outds,
		 Key key, vector<DataPack>& dps)
  {
    assert(dps.size() == (size_t)(NumRegion * NumLattice));

    for (vector<DataPack>::iterator itr = dps.begin(); itr != dps.end();
	 itr++) {
      int newval = *(int*)itr->data->value() + 1;
      Data data(&newval, itr->data->size());
      outds->add(itr->key, data);
    }
    return 0;
  }
};

void run_nicam(DataStore* inds, DataStore* outds)
{
  PseudoNICAM mapper;
  View view(DimEnsembleData);
  bool view_flag[3] = {true, false, false};
  view.set(view_flag);
  inds->map(outds, mapper, view);
}

class PseudoLETKF : public DataStore::Mapper {
public:
  int operator()(DataStore* inds, DataStore* outds,
		 Key key, vector<DataPack>& dps)
  {
    assert(dps.size() == (size_t)NumEnsemble);

    for (vector<DataPack>::iterator itr = dps.begin(); itr != dps.end();
	 itr++) {
      int newval = *(int*)itr->data->value() - 1;
      Data data(&newval, itr->data->size());
      outds->add(itr->key, data);
    }
    return 0;
  }
};

void run_letkf(DataStore* inds, DataStore* outds)
{
  PseudoLETKF mapper;
  View view(DimEnsembleData);
  bool view_flag[3] = {false, true, true};
  view.set(view_flag);
  inds->map(outds, mapper, view);
}
