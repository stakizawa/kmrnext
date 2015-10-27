#include <iostream>
#include "kmrnext.hpp"

using namespace std;

const int Dimension3 = 3;
const int Dim3_0 = 10;
const int Dim3_1 = 10;
const int Dim3_2 = 10;
const int Dimension2 = 2;
const int Dim2_0 = 10;
const int Dim2_1 = 10;

class Loader {
public:
  int operator()(Next::DataStore *ds, const string& file)
  {
    //cout << "From Fanctor: " << file << endl;
    //cout << "From Fanctor: " << ds->to_string() << endl;
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

class Incrementer {
public:
  int operator()(Next::DataStore *inds, Next::DataStore *outds,
		 Next::Key key, vector<Next::DataPack>& dps)
  {
    // cout << "Key: " << key.to_string() << endl;
    // cout << "Data Count: " << dps.size() << endl;
    long sum = 0;
    for (int i = 0; i < dps.size(); i++) {
      Next::DataPack& dp = dps.at(i);
      long v = *(long *)dp.data->value();
      // cout << "  Data: " << dp.key.to_string() << ", " << v << endl;
      sum += v;
    }
    Next::Data d(&sum, sizeof(long));
    outds->add(key, d);
    return 0;
  }
};

class Printer {
public:
  int operator()(Next::DataStore *inds, Next::DataStore *outds,
		 Next::Key key, vector<Next::DataPack>& dps)
  {
    cout << "Key: " << key.to_string() << endl;
    cout << "Data Count: " << dps.size() << endl;
    for (int i = 0; i < dps.size(); i++) {
      Next::DataPack& dp = dps.at(i);
      long v = *(long *)dp.data->value();
      cout << "  Data: " << dp.key.to_string() << ", " << v << endl;
    }
    return 0;
  }
};

void print_gotten_data(vector<Next::DataPack>* dpvec, Next::View& v,
		       Next::Key& k)
{
  cout << "Gotten Data" << endl;
  cout << "  view: " << v.to_string() << endl;
  cout << "  key:  " << k.to_string() << endl;
  cout << "  size: " << dpvec->size() << endl;
  cout << "  values (top10)" << endl;
  int cnt = 0;
  for (vector<Next::DataPack>::iterator itr = dpvec->begin();
       itr != dpvec->end(); itr++) {
    if (cnt >= 10) {
      break;
    }
    cnt += 1;
    cout << "    " << itr->key.to_string() << " : "
	 << *(long *)itr->data->value() << endl;
  }
  cout << endl;
}

int
main()
{
  ///////////  Create a DataStore
  Next::DataStore ds1(Dimension3);
  size_t sizes3[Dimension3] = {Dim3_0, Dim3_1, Dim3_2};
  ds1.set(sizes3);
  cout << "DataStore: " << ds1.to_string() << endl;

  ///////////  Load data contents from a file
  vector<string> files;
  files.push_back("dummy1");
  //  files.push_back("dummy2");
  Loader loader;
  ds1.load_files(files, loader);

  ///////////  Setup keys
  Next::Key key1(Dimension3);
  Next::Key key2(Dimension3);
  size_t kval1[Dimension3] = {2, 2, 2};
  size_t kval2[Dimension3] = {2, 2, 3};
  key1.set(kval1);
  key2.set(kval2);

  ///////////  Get a data from a DataStore
  Next::DataPack dp1 = ds1.get(key1);
  cout << dp1.key.to_string() << " : " << *(long *)dp1.data->value() << endl;
  // cout << "Size: " << dp1.data->size() << endl;
  dp1 = ds1.get(key2);
  cout << dp1.key.to_string() << " : " << *(long *)dp1.data->value() << endl;
  // cout << "Size: " << dp1.data->size() << endl;

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
  vector<Next::DataPack> *dpvec = ds1.get(v1, key1);
  print_gotten_data(dpvec, v1, key1);
  delete dpvec;
  dpvec = ds1.get(v1, key2);
  print_gotten_data(dpvec, v1, key2);
  delete dpvec;
  dpvec = ds1.get(v2, key1);
  print_gotten_data(dpvec, v2, key1);
  delete dpvec;
  dpvec = ds1.get(v2, key2);
  print_gotten_data(dpvec, v2, key2);
  delete dpvec;
  dpvec = ds1.get(v3, key1);
  print_gotten_data(dpvec, v3, key1);
  delete dpvec;
  dpvec = ds1.get(v3, key2);
  print_gotten_data(dpvec, v3, key2);
  delete dpvec;
  dpvec = ds1.get(v4, key1);
  print_gotten_data(dpvec, v4, key1);
  delete dpvec;
  dpvec = ds1.get(v4, key2);
  print_gotten_data(dpvec, v4, key2);
  delete dpvec;

  ///////////  Apply map functions
  Next::DataStore ds2(Dimension2);
  size_t sizes2[Dimension2] = {Dim2_0, Dim2_1};
  ds2.set(sizes2);
  Incrementer incr;
  ds1.map(&ds2, incr, v2);

  ///////////  Get data from ds2
  Next::View v5(Dimension2);
  bool flags5[Dimension2] = {false, false};
  v5.set(flags5);
  Printer printer;
  ds2.map(NULL, printer, v5);

  return 0;
}
