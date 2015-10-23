#include <iostream>
#include <iostream>
#include <sstream>
#include "kmrnext.hpp"

using namespace std;

const int Dimension = 3;
typedef Next::Key<Dimension> Key3;
typedef Next::DataStore<Dimension> DS3;
typedef Next::View<Dimension> V3;

// For convenience
const int Dim0 = 10;
const int Dim1 = 10;
const int Dim2 = 10;

class Loader {
public:
  int operator()(DS3 *ds, const string& file)
  {
    //cout << "From Fanctor: " << file << endl;
    //cout << "From Fanctor: " << ds->to_string() << endl;
    Key3 key;
    for (int i = 0; i < Dim0; i++) {
      key.value[0] = i;
      for (int j = 0; j < Dim1; j++) {
	key.value[1] = j;
	for (int k = 0; k < Dim2; k++) {
	  key.value[2] = k;
	  long val = i*j*k;
	  Next::Data d(&val, sizeof(long));
	  ds->add(key, d);
	}
      }
    }
    return 0;
  }
};

void print_gotten_data(vector<Next::Data>* dvec, V3& v, Key3& k)
{
  cout << "Gotten Data" << endl;
  cout << "  view: " << v.to_string() << endl;
  cout << "  key:  " << k.to_string() << endl;
  cout << "  size: " << dvec->size() << endl;
  cout << "  values (top10)" << endl << "    ";
  int cnt = 0;
  for (vector<Next::Data>::iterator itr = dvec->begin();
       itr != dvec->end(); itr++) {
    if (cnt >= 10) {
      break;
    }
    cnt += 1;
    long val = *(long *)itr->value();
    cout << val << ", ";
  }
  cout << endl;
}

int
main()
{
  cout << "Data dimension: " << DS3::Dimension << endl;

  ///////////  Create a DataStore
  size_t sizes[Dimension] = {Dim0, Dim1, Dim2};
  DS3 ds1(sizes);
  cout << ds1.to_string() << endl;

  ///////////  Load data contents from a file
  vector<string> files;
  files.push_back("dummy1");
  //  files.push_back("dummy2");
  Loader loader;
  ds1.load_files(files, loader);

  ///////////  Setup keys
  size_t kval1[Dimension] = {2, 2, 2};
  size_t kval2[Dimension] = {2, 2, 3};
  Key3 key1(kval1);
  Key3 key2(kval2);

  ///////////  Get a data from a DataStore
  Next::Data d1 = ds1.get(key1);
  cout << "Value: " << *(long *)d1.value() << endl;
  //cout << "Size: " << d1.size() << endl;
  d1 = ds1.get(key2);
  cout << "Value: " << *(long *)d1.value() << endl;
  //cout << "Size: " << d1.size() << endl;

  ///////////  Setup views
  bool flags1[Dimension] = {true, true, true};
  bool flags2[Dimension] = {true, false, true};
  bool flags3[Dimension] = {true, false, false};
  bool flags4[Dimension] = {false, false, false};
  V3 v1(flags1);
  V3 v2(flags2);
  V3 v3(flags3);
  V3 v4(flags4);

  ///////////  Get a data from a DataStore with a view
  vector<Next::Data> *dvec = ds1.get(v1, key1);
  print_gotten_data(dvec, v1, key1);
  delete dvec;
  dvec = ds1.get(v1, key2);
  print_gotten_data(dvec, v1, key2);
  delete dvec;
  dvec = ds1.get(v2, key1);
  print_gotten_data(dvec, v2, key1);
  delete dvec;
  dvec = ds1.get(v2, key2);
  print_gotten_data(dvec, v2, key2);
  delete dvec;
  dvec = ds1.get(v3, key1);
  print_gotten_data(dvec, v3, key1);
  delete dvec;
  dvec = ds1.get(v3, key2);
  print_gotten_data(dvec, v3, key2);
  delete dvec;
  dvec = ds1.get(v4, key1);
  print_gotten_data(dvec, v4, key1);
  delete dvec;
  dvec = ds1.get(v4, key2);
  print_gotten_data(dvec, v4, key2);
  delete dvec;

  ///////////  Apply map functions
  //  ds1.map(v1, mapper);

  return 0;
}
