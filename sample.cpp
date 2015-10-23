#include <iostream>
#include <iostream>
#include <sstream>
#include "kmrnext.hpp"

using namespace std;

const int Dimension = 3;
const int Dim0 = 10;
const int Dim1 = 10;
const int Dim2 = 10;

class Loader {
public:
  int operator()(Next::DataStore *ds, const string& file)
  {
    //cout << "From Fanctor: " << file << endl;
    //cout << "From Fanctor: " << ds->to_string() << endl;
    Next::Key key(Dimension);
    for (int i = 0; i < Dim0; i++) {
      key.set_dim(0, i);
      for (int j = 0; j < Dim1; j++) {
	key.set_dim(1, j);
	for (int k = 0; k < Dim2; k++) {
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
  Next::DataStore ds1(Dimension);
  size_t sizes[Dimension] = {Dim0, Dim1, Dim2};
  ds1.set(sizes);
  cout << "DataStore: " << ds1.to_string() << endl;

  ///////////  Load data contents from a file
  vector<string> files;
  files.push_back("dummy1");
  //  files.push_back("dummy2");
  Loader loader;
  ds1.load_files(files, loader);

  ///////////  Setup keys
  Next::Key key1(Dimension);
  Next::Key key2(Dimension);
  size_t kval1[Dimension] = {2, 2, 2};
  size_t kval2[Dimension] = {2, 2, 3};
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
  Next::View v1(Dimension);
  Next::View v2(Dimension);
  Next::View v3(Dimension);
  Next::View v4(Dimension);
  bool flags1[Dimension] = {true, true, true};
  bool flags2[Dimension] = {true, false, true};
  bool flags3[Dimension] = {true, false, false};
  bool flags4[Dimension] = {false, false, false};
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
  Next::DataStore ds2(Dimension);
  ds1.set(sizes);
  Incrementer incr;
  ds1.map(ds2, incr, v1);

  return 0;
}
