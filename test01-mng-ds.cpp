#include <iostream>
#include "kmrnext.hpp"

using namespace std;

const int DimCell = 2;
const int DimCell0 = 10;
const int DimCell1 = 10;
const size_t DSCellSizes[DimCell] = {DimCell0, DimCell1};

const string File0 = "file0";
const string File1 = "file1";
const string File2 = "file2";
const string File3 = "file3";

void load_file(Next::DataStore& ds, const string& file);

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
  cout << "0. Create DataStores" << endl;
  cout << "  DSCell0: " << dsu0.to_string() << endl;
  cout << "  DSCell1: " << dsu1.to_string() << endl;
  cout << "  DSCell2: " << dsu2.to_string() << endl;
  cout << "  DSCell3: " << dsu3.to_string() << endl;

  /////////// Load data from files
  load_file(dsu0, File0);

  return 0;
}

void load_file(Next::DataStore& ds, const string& file)
{
  vector<string> files;
  files.push_back(file);
  DataLoader loader;
  ds.load_files(files, loader);
}
