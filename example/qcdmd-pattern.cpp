#include <iostream>
#include "kmrnext.hpp"

using namespace std;
using namespace kmrnext;

#if 0
const size_t kEnsembles = 20;
#endif
const size_t kEnsembles = 2;

const size_t kDimSize = 4;
#if 0
const size_t kX = 16;
const size_t kY = 16;
const size_t kZ = 16;
const size_t kT = 16;
#endif
const size_t kX = 2;
const size_t kY = 2;
const size_t kZ = 2;
const size_t kT = 2;

const size_t kQCDSpace[kDimSize] = {kX, kY, kZ, kT};
const size_t kQCDEnsembleSpace[kDimSize + 1] = {kEnsembles, kX, kY, kZ, kT};

const size_t kIterations = 10;

void load_3d_data(DataStore* ds);
void load_4d_data(DataStore* ds);
Key local_key(Key &key, View view);
void print_ds(DataStore* ds);


/// Task class for calculating field C
///
/// It's a dummy class which just increments values in input DS.
class FieldC_Calculator : public DataStore::Mapper {
public:
  int operator()(DataStore* inds, DataStore* outds,
		 Key& key, vector<DataPack>& dps,
		 DataStore::MapEnvironment& env)
  {
    for (vector<DataPack>::iterator itr = dps.begin(); itr != dps.end();
	 itr++) {
      long *data_old = static_cast<long*>(itr->data()->value());
      long data_new = *data_old + 1;
      Data data(&data_new, itr->data()->size());
      outds->add(itr->key(), data);
    }
    return 0;
  }
};

/// Task class for calculating Inv
///
/// It's a dummy class which just increments values in input DS.
class Inv_Calculator : public DataStore::Mapper {
  DataStore* dsC_;
public:
  Inv_Calculator(DataStore *dsC) : dsC_(dsC) {}

  int operator()(DataStore* inds, DataStore* outds,
		 Key& key, vector<DataPack>& dps,
		 DataStore::MapEnvironment& env)
  {
    for (vector<DataPack>::iterator itr = dps.begin(); itr != dps.end();
	 itr++) {
      long *data_old = static_cast<long*>(itr->data()->value());
      Key lkey = local_key(itr->key(), env.view);
      DataPack cdp = dsC_->get(lkey);
      long *cdata = static_cast<long*>(cdp.data()->value());
      long data_new = *data_old + *cdata;
      Data data(&data_new, itr->data()->size());
      outds->add(itr->key(), data);
    }
    return 0;
  }
};

/// Task class for updating P
///
/// It's a dummy class which just increments values in input DS.
class FieldP_Updater : public DataStore::Mapper {
  DataStore* dsC_;
  DataStore* dsG_;
  DataStore* dsF_;
public:
  FieldP_Updater(DataStore *dsC, DataStore *dsG, DataStore *dsF)
    : dsC_(dsC), dsG_(dsG), dsF_(dsF) {}

  int operator()(DataStore* inds, DataStore* outds,
		 Key& key, vector<DataPack>& dps,
		 DataStore::MapEnvironment& env)
  {
    for (vector<DataPack>::iterator itr = dps.begin(); itr != dps.end();
	 itr++) {
      DataPack cdp = dsC_->get(itr->key());
      long *cdata = static_cast<long*>(cdp.data()->value());
      DataPack gdp = dsG_->get(itr->key());
      long *gdata = static_cast<long*>(gdp.data()->value());
      long fdata = 0;
      for (size_t i = 0; i < kEnsembles; i++) {
	Key fkey = Key(kDimSize + 1);
	fkey.set_dim(0, i);
	for (size_t j = 0; j < kDimSize; j++) {
	  fkey.set_dim(j+1, itr->key().dim(j));
	}
	DataPack fdp = dsF_->get(fkey);
	long *fdata0 = static_cast<long*>(fdp.data()->value());
	fdata += *fdata0;
      }

      long *data_old = static_cast<long*>(itr->data()->value());
      long data_new = *data_old - *cdata - *gdata + fdata;
      Data data(&data_new, itr->data()->size());
      outds->add(itr->key(), data);
    }
    return 0;
  }
};

/// Task class for updating G
///
/// It's a dummy class which just increments values in input DS.
class FieldG_Updater : public DataStore::Mapper {
  DataStore* dsP_;
public:
  FieldG_Updater(DataStore *dsP) : dsP_(dsP) {}

  int operator()(DataStore* inds, DataStore* outds,
		 Key& key, vector<DataPack>& dps,
		 DataStore::MapEnvironment& env)
  {
    for (vector<DataPack>::iterator itr = dps.begin(); itr != dps.end();
	 itr++) {
      DataPack pdp = dsP_->get(itr->key());
      long *pdata = static_cast<long*>(pdp.data()->value());

      long *data_old = static_cast<long*>(itr->data()->value());
      long data_new = *data_old + *pdata;
      Data data(&data_new, itr->data()->size());
      outds->add(itr->key(), data);
    }
    return 0;
  }
};


//////////////////////////////////////////////////////////////////////////////
// Main starts from here.
//////////////////////////////////////////////////////////////////////////////
int
main(int argc, char **argv)
{
  KMRNext *next = KMRNext::init(argc, argv);

  // Initialize G
  DataStore* dsG = next->create_ds(kDimSize);
  dsG->set(kQCDSpace);
  load_3d_data(dsG);

  // Initialize P
  DataStore* dsP = next->create_ds(kDimSize);
  dsP->set(kQCDSpace);
  load_3d_data(dsP);

  // Initialize F
  DataStore* dsF0 = next->create_ds(kDimSize + 1);
  dsF0->set(kQCDEnsembleSpace);
  load_4d_data(dsF0);

  print_ds(dsG);
  print_ds(dsP);
  print_ds(dsF0);

  // main loop
  for (size_t i = 0; i < kIterations; i++) {
    // Calculate Field C
    DataStore* dsC = next->create_ds(kDimSize);
    dsC->set(kQCDSpace);
    {
      FieldC_Calculator calculator_c;
      // TODO check if view is appropriate
      View c_calc_view(kDimSize);
      bool flags[kDimSize] = {false, false, false, false};
      c_calc_view.set(flags);
      // TODO set split
      dsG->map(calculator_c, c_calc_view, dsC);
    }

    cout << "Field C" << endl;
    print_ds(dsC);

    // Initialize F1 (output of inv)
    DataStore* dsF1 = next->create_ds(kDimSize + 1);
    dsF1->set(kQCDEnsembleSpace);

    // Calculate inv
    {
      Inv_Calculator calculator_inv(dsC);
      View inv_calc_view(kDimSize + 1);
      bool flags[kDimSize + 1] = {true, false, false, false, false};
      inv_calc_view.set(flags);
      // TODO set split
      dsF0->map(calculator_inv, inv_calc_view, dsF1);
    }

    cout << "Field F'" << endl;
    print_ds(dsF1);

    // Update Field P
    {
      FieldP_Updater updater_p(dsC, dsG, dsF1);
      // TODO check if view is appropriate
      View p_update_view(kDimSize);
      bool flags[kDimSize] = {false, false, false, false};
      p_update_view.set(flags);
      // TODO set split
      dsP->map(updater_p, p_update_view);
    }

    cout << "Field P" << endl;
    print_ds(dsP);

    // Update Field G
    {
      FieldG_Updater updater_g(dsP);
      // TODO check if view is appropriate
      View g_update_view(kDimSize);
      bool flags[kDimSize] = {false, false, false, false};
      g_update_view.set(flags);
      // TODO set split
      dsG->map(updater_g, g_update_view);
    }

    cout << "Field G" << endl;
    print_ds(dsG);

    // Cleanup
    delete dsF1;
    delete dsC;
  }

  delete dsG;
  delete dsP;
  delete dsF0;

  KMRNext::finalize();
  return 0;
}


//////////////////////////////////////////////////////////////////////////////
// Utility functions and classes
//////////////////////////////////////////////////////////////////////////////

class DataLoader3D : public DataStore::Loader<long> {
public:
  int operator()(DataStore* ds, const long& num)
  {
    Key key(kDimSize);
    Data data(static_cast<void*>(const_cast<long*>(&num)), sizeof(long));
    for (size_t i = 0; i < kX; i++) {
      key.set_dim(0, i);
      for (size_t j = 0; j < kY; j++) {
	key.set_dim(1, j);
	for (size_t k = 0; k < kZ; k++) {
	  key.set_dim(2, k);
	  for (size_t l = 0; l < kT; l++) {
	    key.set_dim(3, l);
	    ds->add(key, data);
	  }
	}
      }
    }
    return 0;
  }
};

class DataLoader4D : public DataStore::Loader<long> {
public:
  int operator()(DataStore* ds, const long& num)
  {
    Key key(kDimSize + 1);
    Data data(static_cast<void*>(const_cast<long*>(&num)), sizeof(long));
    for (size_t h = 0; h < kEnsembles; h++) {
      key.set_dim(0, h);
      for (size_t i = 0; i < kX; i++) {
	key.set_dim(1, i);
	for (size_t j = 0; j < kY; j++) {
	  key.set_dim(2, j);
	  for (size_t k = 0; k < kZ; k++) {
	    key.set_dim(3, k);
	    for (size_t l = 0; l < kT; l++) {
	      key.set_dim(4, l);
	      ds->add(key, data);
	    }
	  }
	}
      }
    }
    return 0;
  }
};

void load_data(DataStore* ds, DataStore::Loader<long>& dl) {
  vector<long> data_srcs;
  data_srcs.push_back(0);
  ds->load_integers(data_srcs, dl);
}

void load_3d_data(DataStore* ds) {
  DataLoader3D dl;
  load_data(ds, dl);
}

void load_4d_data(DataStore* ds) {
  DataLoader4D dl;
  load_data(ds, dl);
}

Key local_key(Key &key, View view) {
  size_t ksiz = 0;
  for (size_t i = 0; i < view.size(); i++) {
    if (!view.dim(i)) {
      ksiz += 1;
    }
  }
  Key okey(ksiz);
  size_t kidx = 0;
  for (size_t i = 0; i < view.size(); i++) {
    if (!view.dim(i)) {
      okey.set_dim(kidx, key.dim(i));
      kidx += 1;
    }
  }
  return okey;
}

class DataPrinter : public DataPack::Dumper {
public:
  string operator()(DataPack& dp)
  {
    ostringstream os;
    os << dp.key().to_string() << " : ";
    long *data = static_cast<long*>(dp.data()->value());
    os << *data << endl;
    return os.str();
  }
};

void print_ds(DataStore* ds) {
  DataPrinter dp;
  string str = ds->dump(dp);
  cout << str << endl;
}
