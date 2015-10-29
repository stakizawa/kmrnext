#include <iostream>
#include "kmrnext.hpp"

using namespace std;

const bool print = true;

const int DimCell = 2;
const int DimCell0 = 10;
const int DimCell1 = 10;
const size_t DSCellSizes[DimCell] = {DimCell0, DimCell1};

const string DSUName0 = "DSU0";
const string DSUName1 = "DSU1";
const string DSUName2 = "DSU2";
const string DSUName3 = "DSU3";

const string File0 = "file0";
const string File1 = "file1";
const string File2 = "file2";
const string File3 = "file3";

void load_file(Next::DataStore& ds, const string& file);
void print_ds(Next::DataStore& ds, const Next::View& v, const string& name,
	      int count);


//////////////////////////////////////////////////////////////////////////////
// Main starts from here.
//////////////////////////////////////////////////////////////////////////////
int
main()
{
  /////////// Create DataStores
  Next::DataStore dsu0(DimCell);
  Next::DataStore dsu1(DimCell);
  Next::DataStore dsu2(DimCell);
  Next::DataStore dsu3(DimCell);
  dsu0.set(DSCellSizes);
  dsu1.set(DSCellSizes);
  dsu2.set(DSCellSizes);
  dsu3.set(DSCellSizes);
  if (print) {
    cout << "0. Create DataStores" << endl;
    cout << "  DSCell0: " << dsu0.to_string() << endl;
    cout << "  DSCell1: " << dsu1.to_string() << endl;
    cout << "  DSCell2: " << dsu2.to_string() << endl;
    cout << "  DSCell3: " << dsu3.to_string() << endl;
    cout << endl;
  }

  /////////// Load data from files
  load_file(dsu0, File0);
  load_file(dsu1, File1);
  load_file(dsu2, File2);
  load_file(dsu3, File3);
  if (print) {
    cout << "1. Load data to DataStores" << endl;
    Next::View v(DimCell);
    bool flag[DimCell] = {false, false};
    v.set(flag);
    print_ds(dsu0, v, DSUName0, 10);
    print_ds(dsu1, v, DSUName1, 10);
    print_ds(dsu2, v, DSUName2, 10);
    print_ds(dsu3, v, DSUName3, 10);
    cout << endl;
  }

  /////////// Conmbine DataStores to create a high-dimensinal DataStore
  // 1. Test one DataStore
  Next::DataStore ds0(DimCell + 1);
  vector<Next::DataStore*> dsus;
  dsus.push_back(&dsu0);
  ds0.set_from(dsus);
  if (print) {
    cout << "2. Increase dimension by set_from()" << endl;
    Next::View v(DimCell + 1);
    bool flag[DimCell + 1] = {false, false, false};
    v.set(flag);
    print_ds(ds0, v, "DS_Single", -1);
    cout << endl;
  }
  // 2. Test four DataStores
  Next::DataStore ds1(DimCell + 1);
  dsus.push_back(&dsu1);
  dsus.push_back(&dsu2);
  dsus.push_back(&dsu3);
  ds1.set_from(dsus);
  if (print) {
    cout << "3. Combine four DataStores by set_from()" << endl;
    Next::View v(DimCell + 1);
    bool flag[DimCell + 1] = {false, false, false};
    v.set(flag);
    print_ds(ds1, v, "DS_Four", -1);
    cout << endl;
  }

  /////////// Split a DataStore to low-dimensional DataStores
  Next::DataStore dsm0(DimCell - 1);
  Next::DataStore dsm1(DimCell - 1);
  Next::DataStore dsm2(DimCell - 1);
  Next::DataStore dsm3(DimCell - 1);
  Next::DataStore dsm4(DimCell - 1);
  Next::DataStore dsm5(DimCell - 1);
  Next::DataStore dsm6(DimCell - 1);
  Next::DataStore dsm7(DimCell - 1);
  Next::DataStore dsm8(DimCell - 1);
  Next::DataStore dsm9(DimCell - 1);
  vector<Next::DataStore*> dsms;
  dsms.push_back(&dsm0);
  dsms.push_back(&dsm1);
  dsms.push_back(&dsm2);
  dsms.push_back(&dsm3);
  dsms.push_back(&dsm4);
  dsms.push_back(&dsm5);
  dsms.push_back(&dsm6);
  dsms.push_back(&dsm7);
  dsms.push_back(&dsm8);
  dsms.push_back(&dsm9);
  dsu0.split_to(dsms);
  if (print) {
    cout << "4. Split a DataStore to 10 DataStores by split_to()"  << endl;
    Next::View v(DimCell - 1);
    bool flag[DimCell - 1] = {false};
    v.set(flag);
    print_ds(*dsms.at(0), v, "DS_Smallest", -1);
    print_ds(*dsms.at(1), v, "DS_Smallest", -1);
    cout << endl;
  }

  return 0;
}


class DataLoader {
public:
  int operator()(Next::DataStore *ds, const string& file)
  {
    Next::Key key(DimCell);
    for (int i = 0; i < DimCell0; i++) {
      key.set_dim(0, i);
      for (int j = 0; j < DimCell1; j++) {
    	key.set_dim(1, j);
	long val;
	if (file == File0) {
	  val = 1;
	} else if (file == File1) {
	  val = 2;
	} else if (file == File2) {
	  val = 3;
	} else if (file == File3) {
	  val = 4;
	} else {
	  val = 0;
	}
	Next::Data d(&val, sizeof(long));
	ds->add(key, d);
      }
    }
    return 0;
  }
};

// A mapper class that prints all data.
class DataStorePrinter {
  int _max_count;
  string _padding;

public:
  DataStorePrinter(const int max_count, const string& padding)
    : _max_count(max_count), _padding(padding) {}

  int operator()(Next::DataStore *inds, Next::DataStore *outds,
		 Next::Key key, vector<Next::DataPack>& dps)
  {
    //cout << _padding << "Key: " << key.to_string() << endl;
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

void load_file(Next::DataStore& ds, const string& file)
{
  vector<string> files;
  files.push_back(file);
  DataLoader loader;
  ds.load_files(files, loader);
}

void print_ds(Next::DataStore& ds, const Next::View& v, const string& name,
	      int count)
{
  cout << "  " << name << endl;
  DataStorePrinter dsp(count, "    ");
  ds.map(NULL, dsp, v);
}
