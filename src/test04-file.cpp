#include <iostream>
#include "kmrnext.hpp"
#ifdef BACKEND_KMR
#include <mpi.h>
#endif
/// This is just a simplified version of test00-basic.cpp and
/// test01-mng-ds.cpp, except that Data in DataStore are stored
/// in files.

using namespace std;

int rank = 0;

const bool kPrint = true;

const size_t kDim3 = 3;
const size_t kDim3_0 = 2;
const size_t kDim3_1 = 4;
const size_t kDim3_2 = 4;
const size_t kDim2 = 2;
const size_t kDim2_0 = 4;
const size_t kDim2_1 = 4;
const size_t kDim1 = 1;
const size_t kDim1_0 = 4;


void load_data(kmrnext::DataStore* ds);
void print_line(string str);
void print_ds_count(kmrnext::DataStore* ds, string padding="    ");
void print_ds_contents(kmrnext::DataStore* ds, string padding="    ");
void print_datapack(kmrnext::DataPack& dp, string padding="    ");
void print_datapacks(vector<kmrnext::DataPack>* dpvec, string padding="    ");

class Summarizer : public kmrnext::DataStore::Mapper {
public:
  int operator()(kmrnext::DataStore* inds, kmrnext::DataStore* outds,
		 kmrnext::Key& key, vector<kmrnext::DataPack>& dps,
		 kmrnext::DataStore::MapEnvironment& env) {
    long sum = 0;
    for (size_t i = 0; i < dps.size(); i++) {
      kmrnext::DataPack& dp = dps.at(i);
      long v = *static_cast<long*>(dp.data().value());
      sum += v;
    }
    kmrnext::Data d(&sum, sizeof(long));
    outds->add(key, d);
    return 0;
  }
};

class OneSetter : public kmrnext::DataStore::Mapper {
public:
  int operator()(kmrnext::DataStore* inds, kmrnext::DataStore* outds,
		 kmrnext::Key& key, vector<kmrnext::DataPack>& dps,
		 kmrnext::DataStore::MapEnvironment& env) {
    long val = 1;
    kmrnext::Data d(&val, sizeof(long));
    outds->add(key, d);
    return 0;
  }
};


//////////////////////////////////////////////////////////////////////////////
// Main starts from here.
//////////////////////////////////////////////////////////////////////////////
int
main(int argc, char **argv) {
  kmrnext::KMRNext *next = kmrnext::KMRNext::init(argc, argv);
  next->set_io_mode(kmrnext::KMRNext::File);

#ifdef BACKEND_KMR
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  ///////////  Create a DataStore
  kmrnext::DataStore *ds0 = next->create_ds(kDim2);
  size_t size2[kDim2] = {kDim2_0, kDim2_1};
  ds0->set(size2);
  if (kPrint) {
    print_line("0. Just creating a DataStore");
    print_line("    DataStore: " + ds0->to_string());
    print_ds_count(ds0, "        ");
  }

  ///////////  Load data contents
  load_data(ds0);
  if (kPrint) {
    print_line("1. Load Data to the DataStore");
    print_ds_count(ds0);
    print_ds_contents(ds0, "        ");
    print_line("");
  }

  ///////////  Get a data from the DataStore
  kmrnext::Key key0(kDim2);
  size_t kval0[kDim2] = {2, 1};
  key0.set(kval0);
  kmrnext::DataPack dp0 = ds0->get(key0);
  if (kPrint) {
    print_line("2. Get a Data from the DataStore");
    print_line("    Key:  " + key0.to_string());
    print_datapack(dp0, "        ");
  }
  kmrnext::Key key1(kDim2);
  size_t kval1[kDim2] = {3, 3};
  key1.set(kval1);
  kmrnext::DataPack dp1 = ds0->get(key1);
  if (kPrint) {
    print_line("    Key:  " + key1.to_string());
    print_datapack(dp1, "        ");
    print_line("");
  }

  ///////////  Get multiple data from the DataStore
  kmrnext::View v0(kDim2);
  bool v0flag[kDim2] = {true, false};
  v0.set(v0flag);
  vector<kmrnext::DataPack> *dpvec0 = ds0->get(v0, key0);
  if (kPrint) {
    print_line("3. Get multiple Data from the DataStore");
    print_line("    Key:  " + key0.to_string());
    print_line("    View: " + v0.to_string());
    print_datapacks(dpvec0, "        ");
  }
  delete dpvec0;
  kmrnext::View v1(kDim2);
  bool v1flag[kDim2] = {false, true};
  v1.set(v1flag);
  vector<kmrnext::DataPack> *dpvec1 = ds0->get(v1, key1);
  if (kPrint) {
    print_line("    Key:  " + key1.to_string());
    print_line("    View: " + v1.to_string());
    print_datapacks(dpvec1, "        ");
    print_line("");
  }
  delete dpvec1;

  ///////////  Remove and then add a Data
  kmrnext::DataPack dp2 = ds0->remove(key0);
  if (kPrint) {
    print_line("4. Remove a Data and add a new one");
    print_line("    Remove:  " + key0.to_string());
    print_datapack(dp2, "        ");
    print_line("    DataStore: " + ds0->to_string());
    print_ds_count(ds0, "        ");
    print_ds_contents(ds0, "        ");
    print_line("");
  }
  kmrnext::DataPack dp3 = ds0->remove(key1);
  if (kPrint) {
    print_line("    Remove:  " + key1.to_string());
    print_datapack(dp3, "        ");
    print_line("    DataStore: " + ds0->to_string());
    print_ds_count(ds0, "        ");
    print_ds_contents(ds0, "        ");
    print_line("");
  }
  ds0->add(key0, dp3.data());
  ds0->add(key1, dp2.data());
  if (kPrint) {
    print_line("    Add:  " + key0.to_string());
    print_line("    Add:  " + key1.to_string());
    print_line("    DataStore: " + ds0->to_string());
    print_ds_count(ds0, "        ");
    print_ds_contents(ds0, "        ");
    print_line("");
  }

  ///////////  Duplicate the DataStore
  kmrnext::DataStore *ds1 = ds0->duplicate();
  if (kPrint) {
    print_line("5. Duplicate the DataStore");
    print_line("    Duplicated DataStore: " + ds1->to_string());
    print_ds_count(ds1);
    print_ds_contents(ds1, "        ");
    print_line("");
  }

  ///////////  Map the DataStore
  kmrnext::DataStore *ds2 = next->create_ds(kDim1);
  size_t size1[kDim1] = {kDim1_0};
  ds2->set(size1);
  Summarizer sumr;
  ds0->map(sumr, v0, ds2);
  if (kPrint) {
    print_line("6. Map the DataStore");
    print_line("    DataStore: " + ds0->to_string());
    print_ds_count(ds2);
    print_ds_contents(ds2, "        ");
    print_line("");
  }
  delete ds2;

  ///////////  Map the DataStore in-place
  kmrnext::View v2(kDim2);
  bool v2flag[kDim2] = {true, true};
  v2.set(v2flag);
  OneSetter oner;
  ds1->map(oner, v2);
  if (kPrint) {
    print_line("7. Map the DataStore in-place");
    print_line("    DataStore: " + ds1->to_string());
    print_ds_count(ds1);
    print_ds_contents(ds1, "        ");
    print_line("");
  }

  ///////////  Merge two DataStores
  kmrnext::DataStore *ds3 = next->create_ds(kDim3);
  vector<kmrnext::DataStore*> indslst;
  indslst.push_back(ds0);
  indslst.push_back(ds1);
  ds3->set_from(indslst);
  if (kPrint) {
    print_line("8. Merge two DataStores");
    print_line("    DataStore: " + ds3->to_string());
    print_ds_count(ds3);
    print_ds_contents(ds3, "        ");
    print_line("");
  }
  kmrnext::Key key2(kDim3);
  size_t kval2[kDim3] = {0, 1, 2};
  key2.set(kval2);
  kmrnext::DataPack dp4 = ds3->get(key2);
  if (kPrint) {
    print_line("    Key:  " + key2.to_string());
    print_datapack(dp4, "        ");
    print_line("");
  }
  delete ds3;
  delete ds1;

  ///////////  Split into four DataStores
  kmrnext::DataStore *sds0 = next->create_ds(kDim1);
  kmrnext::DataStore *sds1 = next->create_ds(kDim1);
  kmrnext::DataStore *sds2 = next->create_ds(kDim1);
  kmrnext::DataStore *sds3 = next->create_ds(kDim1);
  vector<kmrnext::DataStore*> outdslst;
  outdslst.push_back(sds0);
  outdslst.push_back(sds1);
  outdslst.push_back(sds2);
  outdslst.push_back(sds3);
  ds0->split_to(outdslst);
  if (kPrint) {
    print_line("8. Split into four DataStores");
    print_line("    DataStore 0: " + sds0->to_string());
    print_ds_count(sds0);
    print_ds_contents(sds0, "        ");
    print_line("");
    print_line("    DataStore 1: " + sds1->to_string());
    print_ds_count(sds1);
    print_ds_contents(sds1, "        ");
    print_line("");
    print_line("    DataStore 2: " + sds2->to_string());
    print_ds_count(sds2);
    print_ds_contents(sds2, "        ");
    print_line("");
    print_line("    DataStore 3: " + sds3->to_string());
    print_ds_count(sds3);
    print_ds_contents(sds3, "        ");
    print_line("");
  }
  delete sds0;
  delete sds1;
  delete sds2;
  delete sds3;

  delete ds0;
  kmrnext::KMRNext::finalize();
  return 0;
}


class IntLoader2 : public kmrnext::DataStore::Loader<long> {
public:
  int operator()(kmrnext::DataStore* ds, const long& i) {
    kmrnext::Key key(kDim2);
    key.set_dim(0, i);
    for (size_t j = 0; j < kDim2_1; j++) {
      key.set_dim(1, j);
      long val = static_cast<long>(i * kDim2_1 + j);
      kmrnext::Data d(&val, sizeof(long));
      ds->add(key, d);
    }
    return 0;
  }
};

void load_data(kmrnext::DataStore* ds) {
  vector<long> ints;
  ints.push_back(0);
  ints.push_back(1);
  ints.push_back(2);
  ints.push_back(3);
  IntLoader2 loader;
  ds->load_integers(ints, loader);
}

void print_line(string str) {
  if (rank != 0) return;
  cout << str << endl;
}

void print_ds_count(kmrnext::DataStore* ds, string padding) {
  long count = ds->count();
  ostringstream os;
  os << padding << "Count: " << count;
  print_line(os.str());
}

class DPPrinter : public kmrnext::DataPack::Dumper {
public:
  string operator()(kmrnext::DataPack& dp) {
    ostringstream os;
    os << dp.key().to_string() << " : "
       << *static_cast<long*>(dp.data().value()) << endl;
    return os.str();
  }
};

void print_ds_contents(kmrnext::DataStore* ds, string padding) {
  DPPrinter printer;
  string ds_str = ds->dump(printer);
  if (rank != 0) return;
  istringstream is(ds_str);
  string str;
  while(getline(is, str, '\n')) {
    cout << padding << str << endl;
  }
}

void print_datapack(kmrnext::DataPack& dp, string padding) {
  if (rank != 0) return;
  kmrnext::Key key = dp.key();
  kmrnext::Data data = dp.data();
  cout << padding << key.to_string() << " : "
       << *static_cast<long*>(data.value()) << endl;
}

void print_datapacks(vector<kmrnext::DataPack>* dpvec, string padding) {
  if (rank != 0) return;
  for (vector<kmrnext::DataPack>::iterator itr = dpvec->begin();
       itr != dpvec->end(); itr++) {
    kmrnext::Key key = itr->key();
    kmrnext::Data data = itr->data();
    cout << padding << key.to_string() << " : "
	 << *static_cast<long*>(data.value()) << endl;
  }
}
