#include <iostream>
#include "kmrnext.hpp"

using namespace std;

const bool print = true;

const int Dimension3 = 3;
const int Dim3_0 = 10;
const int Dim3_1 = 10;
const int Dim3_2 = 10;
const int Dimension2 = 2;
const int Dim2_0 = 10;
const int Dim2_1 = 10;

void load_data(Next::DataStore& ds);
void print_get_result(Next::Key& key, Next::DataPack& dp);
void print_get_view_result(vector<Next::DataPack>* dpvec, Next::View& v,
			   Next::Key& k, int count);

// A mapper class that calculates sum of data.
class Summarizer : public Next::DataStore::Mapper {
public:
  int operator()(Next::DataStore *inds, Next::DataStore *outds,
		 Next::Key key, vector<Next::DataPack>& dps)
  {
    long sum = 0;
    for (size_t i = 0; i < dps.size(); i++) {
      Next::DataPack& dp = dps.at(i);
      long v = *(long *)dp.data->value();
      sum += v;
    }
    Next::Data d(&sum, sizeof(long));
    outds->add(key, d);
    return 0;
  }
};

// A mapper class that prints all data.
class DataStorePrinter : public Next::DataStore::Mapper {
  int _max_count;
  string _padding;

public:
  DataStorePrinter(const int max_count, const string& padding)
    : _max_count(max_count), _padding(padding) {}

  int operator()(Next::DataStore *inds, Next::DataStore *outds,
		 Next::Key key, vector<Next::DataPack>& dps)
  {
    cout << _padding << "Key: " << key.to_string() << endl;
    cout << _padding << "Count: " << dps.size() << endl;
    if (_max_count > 0) {
      cout << _padding << "Values (top" << _max_count << ")" << endl;
    } else {
      cout << _padding << "Values (all)" << endl;
    }
    int count = 0;
    for (vector<Next::DataPack>::iterator itr = dps.begin();
    	 itr != dps.end(); itr++) {
      if (_max_count > 0 && count++ >= _max_count) {
    	break;
      }
      cout << _padding << "  " << itr->key.to_string() << " : "
    	   << *(long *)itr->data->value() << endl;
    }
    return 0;
  }
};


//////////////////////////////////////////////////////////////////////////////
// Main starts from here.
//////////////////////////////////////////////////////////////////////////////
int
main()
{
  ///////////  Create a DataStore
  Next::DataStore ds1(Dimension3);
  size_t sizes3[Dimension3] = {Dim3_0, Dim3_1, Dim3_2};
  ds1.set(sizes3);
  cout << "DataStore: " << ds1.to_string() << endl;

  ///////////  Load data contents from a file
  load_data(ds1);

  ///////////  Setup keys
  Next::Key key1(Dimension3);
  Next::Key key2(Dimension3);
  size_t kval1[Dimension3] = {2, 2, 2};
  size_t kval2[Dimension3] = {2, 2, 3};
  key1.set(kval1);
  key2.set(kval2);

  ///////////  Get a data from a DataStore
  Next::DataPack dp1 = ds1.get(key1);
  Next::DataPack dp2 = ds1.get(key2);
  if (print) {
    cout << "1. Get a data from a DataStore by get()" << endl;
    print_get_result(key1, dp1);
    print_get_result(key2, dp2);
    cout << endl;
  }

  ///////////  Setup views
  Next::View v1(Dimension3);
  Next::View v2(Dimension3);
  Next::View v3(Dimension3);
  Next::View v4(Dimension3);
  bool flags1[Dimension3] = {true, true, true};
  bool flags2[Dimension3] = {true, false, true};
  bool flags3[Dimension3] = {true, false, false};
  bool flags4[Dimension3] = {false, false, false};
  v1.set(flags1);
  v2.set(flags2);
  v3.set(flags3);
  v4.set(flags4);

  ///////////  Get a data from a DataStore with a view
  vector<Next::DataPack> *dpvec1 = ds1.get(v1, key1);
  vector<Next::DataPack> *dpvec2 = ds1.get(v1, key2);
  vector<Next::DataPack> *dpvec3 = ds1.get(v2, key1);
  vector<Next::DataPack> *dpvec4 = ds1.get(v2, key2);
  vector<Next::DataPack> *dpvec5 = ds1.get(v3, key1);
  vector<Next::DataPack> *dpvec6 = ds1.get(v3, key2);
  vector<Next::DataPack> *dpvec7 = ds1.get(v4, key1);
  vector<Next::DataPack> *dpvec8 = ds1.get(v4, key2);
  if (print) {
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
  Next::DataStore ds2(Dimension2);
  size_t sizes2[Dimension2] = {Dim2_0, Dim2_1};
  ds2.set(sizes2);
  Summarizer sumr;
  ds1.map(&ds2, sumr, v2);
  if (print) {
    cout << "3. Apply map to each data in a DataStore" << endl;
    Next::View v5(Dimension2);
    bool flags5[Dimension2] = {false, false};
    v5.set(flags5);
    DataStorePrinter printer(-1, "  ");
    ds2.map(NULL, printer, v5);
    cout << endl;
  }

  return 0;
}


class Loader : public Next::DataStore::Loader<string> {
public:
  int operator()(Next::DataStore *ds, const string& file)
  {
    Next::Key key(Dimension3);
    for (int i = 0; i < Dim3_0; i++) {
      key.set_dim(0, i);
      for (int j = 0; j < Dim3_1; j++) {
	key.set_dim(1, j);
	for (int k = 0; k < Dim3_2; k++) {
	  key.set_dim(2, k);
	  long val = i*j*k;
	  Next::Data d(&val, sizeof(long));
	  ds->add(key, d);
	}
      }
    }
    return 0;
  }
};

void load_data(Next::DataStore& ds)
{
  vector<string> files;
  files.push_back("dummy1");
  //  files.push_back("dummy2");
  Loader loader;
  ds.load_files(files, loader);
}

void print_get_result(Next::Key& key, Next::DataPack& dp)
{
  cout << "  Query key: " << key.to_string()
       << "    Result: " << dp.key.to_string() << ": "
       << *(long *)dp.data->value()
       << " (Size:" << dp.data->size() << ")" << endl;
}

void print_get_view_result(vector<Next::DataPack>* dpvec, Next::View& v,
			   Next::Key& k, int count)
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
  for (vector<Next::DataPack>::iterator itr = dpvec->begin();
       itr != dpvec->end(); itr++) {
    if (count > 0 && cnt >= count) {
      break;
    }
    cnt += 1;
    cout << "    " << itr->key.to_string() << " : "
  	 << *(long *)itr->data->value() << endl;
  }
  cout << endl;
}
