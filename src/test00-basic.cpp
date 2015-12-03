#include <iostream>
#include "kmrnext.hpp"
#ifdef BACKEND_KMR
#include <mpi.h>
#endif

using namespace std;

int rank = 0;

const bool kPrint = true;

const size_t kDimension3 = 3;
const size_t kDim3_0 = 10;
const size_t kDim3_1 = 10;
const size_t kDim3_2 = 10;
const size_t kDimension2 = 2;
const size_t kDim2_0 = 10;
const size_t kDim2_1 = 10;

void load_data(kmrnext::DataStore *ds);
bool prints();
void print_data_store(kmrnext::DataStore *ds);
void print_get_result(kmrnext::Key& key, kmrnext::DataPack& dp);
void print_get_view_result(vector<kmrnext::DataPack>* dpvec, kmrnext::View& v,
			   kmrnext::Key& k, int count);

// A mapper class that calculates sum of data.
class Summarizer : public kmrnext::DataStore::Mapper {
public:
  int operator()(kmrnext::DataStore *inds, kmrnext::DataStore *outds,
		 kmrnext::Key& key, vector<kmrnext::DataPack>& dps,
		 kmrnext::DataStore::MapEnvironment& env)
  {
    long sum = 0;
    for (size_t i = 0; i < dps.size(); i++) {
      kmrnext::DataPack& dp = dps.at(i);
      long v = *(long *)dp.data()->value();
      sum += v;
    }
    kmrnext::Data d(&sum, sizeof(long));
    outds->add(key, d);
    return 0;
  }
};

// A mapper class that prints all data.
class DataStorePrinter : public kmrnext::DataStore::Mapper {
  int _max_count;
  string _padding;

public:
  DataStorePrinter(const int max_count, const string& padding)
    : _max_count(max_count), _padding(padding) {}

  int operator()(kmrnext::DataStore *inds, kmrnext::DataStore *outds,
		 kmrnext::Key& key, vector<kmrnext::DataPack>& dps,
		 kmrnext::DataStore::MapEnvironment& env)
  {
    cout << _padding << "Key: " << key.to_string() << endl;
    cout << _padding << "Count: " << dps.size() << endl;
    if (_max_count > 0) {
      cout << _padding << "Values (top" << _max_count << ")" << endl;
    } else {
      cout << _padding << "Values (all)" << endl;
    }
    int count = 0;
    for (vector<kmrnext::DataPack>::iterator itr = dps.begin();
    	 itr != dps.end(); itr++) {
      if (_max_count > 0 && count++ >= _max_count) {
    	break;
      }
      cout << _padding << "  " << itr->key().to_string() << " : "
    	   << *(long *)itr->data()->value() << endl;
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
  kmrnext::KMRNext *next = kmrnext::KMRNext::init(argc, argv);

#ifdef BACKEND_KMR
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  ///////////  Create a DataStore
  kmrnext::DataStore *ds1 = next->create_ds(kDimension3);
  size_t sizes3[kDimension3] = {kDim3_0, kDim3_1, kDim3_2};
  ds1->set(sizes3);
  print_data_store(ds1);

  ///////////  Load data contents from a file
  load_data(ds1);

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
  if (prints()) {
    cout << "1. Get a data from a DataStore by get()" << endl;
    print_get_result(key1, dp1);
    print_get_result(key2, dp2);
    cout << endl;
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
  if (prints()) {
    cout << "2. Get data from a DataStore by get(view)" << endl;
    print_get_view_result(dpvec1, v1, key1, 10);
    print_get_view_result(dpvec2, v1, key2, 10);
    print_get_view_result(dpvec3, v2, key1, 10);
    print_get_view_result(dpvec4, v2, key2, 10);
    print_get_view_result(dpvec5, v3, key1, 10);
    print_get_view_result(dpvec6, v3, key2, 10);
    print_get_view_result(dpvec7, v4, key1, 10);
    print_get_view_result(dpvec8, v4, key2, 10);
    cout << endl;
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
  ds1->map(ds2, sumr, v2);
  if (prints()) {
    cout << "3. Apply map to each data in a DataStore" << endl;
    kmrnext::View v5(kDimension2);
    bool flags5[kDimension2] = {false, false};
    v5.set(flags5);
    DataStorePrinter printer(-1, "  ");
    ds2->map(NULL, printer, v5);
    cout << endl;
  }
  delete ds2;

  delete ds1;
  kmrnext::KMRNext::finalize();
  return 0;
}


class Loader : public kmrnext::DataStore::Loader<string> {
public:
  int operator()(kmrnext::DataStore *ds, const string& file)
  {
    kmrnext::Key key(kDimension3);
    for (size_t i = 0; i < kDim3_0; i++) {
      key.set_dim(0, i);
      for (size_t j = 0; j < kDim3_1; j++) {
	key.set_dim(1, j);
	for (size_t k = 0; k < kDim3_2; k++) {
	  key.set_dim(2, k);
	  long val = (long)(i*j*k);
	  kmrnext::Data d(&val, sizeof(long));
	  ds->add(key, d);
	}
      }
    }
    return 0;
  }
};

void load_data(kmrnext::DataStore *ds)
{
  vector<string> files;
  files.push_back("dummy1");
  //  files.push_back("dummy2");
  Loader loader;
  ds->load_files(files, loader);
}

bool prints()
{
  return (rank == 0) && kPrint;
}

void print_data_store(kmrnext::DataStore *ds)
{
  if (rank == 0) {
    cout << "DataStore: " << ds->to_string() << endl;
  }
}

void print_get_result(kmrnext::Key& key, kmrnext::DataPack& dp)
{
  cout << "  Query key: " << key.to_string()
       << "    Result: " << dp.key().to_string() << ": "
       << *(long *)dp.data()->value()
       << " (Size:" << dp.data()->size() << ")" << endl;
}

void print_get_view_result(vector<kmrnext::DataPack>* dpvec, kmrnext::View& v,
			   kmrnext::Key& k, int count)
{
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
  	 << *(long *)itr->data()->value() << endl;
  }
  cout << endl;
}
