#include <iostream>
#include "kmrnext.hpp"
#ifdef BACKEND_KMR
#include <mpi.h>
#endif

using namespace std;

int rank = 0;

const bool kPrint = true;
const int  kDumpCount = 20;

const size_t kDimension3 = 3;
const size_t kDim3_0 = 10;
const size_t kDim3_1 = 10;
const size_t kDim3_2 = 10;
const size_t kDimension2 = 2;
const size_t kDim2_0 = 10;
const size_t kDim2_1 = 10;

void load_data(kmrnext::DataStore* ds);
void print_line(string str);
void print_data_store(kmrnext::DataStore* ds, string padding, int count);
void print_get_result(kmrnext::Key& key, kmrnext::DataPack& dp);
void print_get_view_result(vector<kmrnext::DataPack>* dpvec, kmrnext::View& v,
			   kmrnext::Key& k, int count);

// A mapper class that calculates sum of data.
class Summarizer : public kmrnext::DataStore::Mapper {
public:
  int operator()(kmrnext::DataStore* inds, kmrnext::DataStore* outds,
		 kmrnext::Key& key, vector<kmrnext::DataPack>& dps,
		 kmrnext::DataStore::MapEnvironment& env) {
    long sum = 0;
    for (size_t i = 0; i < dps.size(); i++) {
      kmrnext::DataPack& dp = dps.at(i);
      long v = *static_cast<long*>(dp.data()->value());
      sum += v;
    }
    kmrnext::Data d(&sum, sizeof(long));
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

#ifdef BACKEND_KMR
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  ///////////  Create a DataStore
  kmrnext::DataStore *ds1 = next->create_ds(kDimension3);
  size_t sizes3[kDimension3] = {kDim3_0, kDim3_1, kDim3_2};
  ds1->set(sizes3);
  if (kPrint) {
    print_line("0. Create a DataStore");
    print_line("  DataStore: " + ds1->to_string());
    print_line("");
  }

  ///////////  Load data contents from a file
  load_data(ds1);
  if (kPrint) {
    print_line("1. Load data to a DataStore");
    print_data_store(ds1, "  ", kDumpCount);
    print_line("");
  }

  ///////////  Setup keys
  kmrnext::Key key1(kDimension3);
  kmrnext::Key key2(kDimension3);
  size_t kval1[kDimension3] = {2, 2, 2};
  size_t kval2[kDimension3] = {2, 2, 3};
  key1.set(kval1);
  key2.set(kval2);

  ///////////  Get a data from a DataStore
  kmrnext::DataPack dp1 = ds1->get(key1);
  kmrnext::DataPack dp2 = ds1->get(key2);
  if (kPrint) {
    print_line("2. Get a data from a DataStore by get()");
    print_get_result(key1, dp1);
    print_get_result(key2, dp2);
    print_line("");
  }

  ///////////  Setup views
  kmrnext::View v1(kDimension3);
  kmrnext::View v2(kDimension3);
  kmrnext::View v3(kDimension3);
  kmrnext::View v4(kDimension3);
  bool flags1[kDimension3] = {true, true, true};
  bool flags2[kDimension3] = {true, false, true};
  bool flags3[kDimension3] = {true, false, false};
  bool flags4[kDimension3] = {false, false, false};
  v1.set(flags1);
  v2.set(flags2);
  v3.set(flags3);
  v4.set(flags4);

  ///////////  Get a data from a DataStore with a view
  vector<kmrnext::DataPack> *dpvec1 = ds1->get(v1, key1);
  vector<kmrnext::DataPack> *dpvec2 = ds1->get(v1, key2);
  vector<kmrnext::DataPack> *dpvec3 = ds1->get(v2, key1);
  vector<kmrnext::DataPack> *dpvec4 = ds1->get(v2, key2);
  vector<kmrnext::DataPack> *dpvec5 = ds1->get(v3, key1);
  vector<kmrnext::DataPack> *dpvec6 = ds1->get(v3, key2);
  vector<kmrnext::DataPack> *dpvec7 = ds1->get(v4, key1);
  vector<kmrnext::DataPack> *dpvec8 = ds1->get(v4, key2);
  if (kPrint) {
    print_line("3. Get data from a DataStore by get(view)");
    print_get_view_result(dpvec1, v1, key1, kDumpCount);
    print_get_view_result(dpvec2, v1, key2, kDumpCount);
    print_get_view_result(dpvec3, v2, key1, kDumpCount);
    print_get_view_result(dpvec4, v2, key2, kDumpCount);
    print_get_view_result(dpvec5, v3, key1, kDumpCount);
    print_get_view_result(dpvec6, v3, key2, kDumpCount);
    print_get_view_result(dpvec7, v4, key1, kDumpCount);
    print_get_view_result(dpvec8, v4, key2, kDumpCount);
    print_line("");
  }
  delete dpvec1;
  delete dpvec2;
  delete dpvec3;
  delete dpvec4;
  delete dpvec5;
  delete dpvec6;
  delete dpvec7;
  delete dpvec8;

  ///////////  Apply map functions
  kmrnext::DataStore *ds2 = next->create_ds(kDimension2);
  size_t sizes2[kDimension2] = {kDim2_0, kDim2_1};
  ds2->set(sizes2);
  Summarizer sumr;
  ds1->map(sumr, v2, ds2);
  if (kPrint) {
    print_line("4. Apply map to each data in a DataStore");
    print_line("  Output DataStore");
    print_data_store(ds2, "    ", kDumpCount);
    print_line("");
  }
  delete ds2;

  delete ds1;
  kmrnext::KMRNext::finalize();
  return 0;
}


class IntLoader : public kmrnext::DataStore::Loader<long> {
public:
  int operator()(kmrnext::DataStore* ds, const long& i) {
    kmrnext::Key key(kDimension3);
    key.set_dim(0, i);
    for (size_t j = 0; j < kDim3_1; j++) {
      key.set_dim(1, j);
      for (size_t k = 0; k < kDim3_2; k++) {
	key.set_dim(2, k);
	long val = static_cast<long>(i*j*k);
	kmrnext::Data d(&val, sizeof(long));
	ds->add(key, d);
      }
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
  ints.push_back(4);
  ints.push_back(5);
  ints.push_back(6);
  ints.push_back(7);
  ints.push_back(8);
  ints.push_back(9);
  IntLoader loader;
  ds->load_integers(ints, loader);
}

void print_line(string str) {
  if (rank != 0) return;
  cout << str << endl;
}

// A class for dumping a DataPack
class DPPrinter : public kmrnext::DataPack::Dumper {
public:
  string operator()(kmrnext::DataPack& dp) {
    ostringstream os;
    os << "  " << dp.key().to_string() << " : "
       << *static_cast<long*>(dp.data()->value()) << endl;
    return os.str();
  }
};

void print_data_store(kmrnext::DataStore* ds, string padding, int count) {
  long ds_count = ds->count();
  DPPrinter printer;
  string ds_str = ds->dump(printer);
  if (rank != 0) return;

  cout << padding << "Count of data in the DataStore: " << ds_count << endl;
  if (count > 0) {
    cout << padding << "Values (top" << count << ")" << endl;
  } else {
    cout << padding << "Values (all)" << endl;
  }
  int cnt = 0;
  istringstream is(ds_str);
  string str;
  while(getline(is, str, '\n')) {
    if (count > 0 && cnt >= count) {
      break;
    }
    cnt += 1;
    cout << padding << str << endl;
  }
}

void print_get_result(kmrnext::Key& key, kmrnext::DataPack& dp) {
  if (rank != 0) return;
  cout << "  Query key: " << key.to_string()
       << "    Result: " << dp.key().to_string() << ": "
       << *static_cast<long*>(dp.data()->value())
       << " (Size:" << dp.data()->size() << ")" << endl;
}

void print_get_view_result(vector<kmrnext::DataPack>* dpvec, kmrnext::View& v,
			   kmrnext::Key& k, int count) {
  if (rank != 0) return;
  cout << "  Condition" << endl;
  cout << "    view: " << v.to_string() << ", key: " << k.to_string() << endl;
  cout << "  Result" << endl;
  cout << "    size: " << dpvec->size() << endl;
  if (count > 0) {
    cout << "    values (top" << count << ")" << endl;
  } else {
    cout << "    values (all)" << endl;
  }
  int cnt = 0;
  for (vector<kmrnext::DataPack>::iterator itr = dpvec->begin();
       itr != dpvec->end(); itr++) {
    if (count > 0 && cnt >= count) {
      break;
    }
    cnt += 1;
    cout << "    " << itr->key().to_string() << " : "
  	 << *static_cast<long*>(itr->data()->value()) << endl;
  }
  cout << endl;
}
