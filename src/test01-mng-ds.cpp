#include <iostream>
#include "kmrnext.hpp"

using namespace std;

const bool kPrint = true;

const size_t kDimCell = 2;
const size_t kDimCell0 = 10;
const size_t kDimCell1 = 10;
const size_t kDSCellSizes[kDimCell] = {kDimCell0, kDimCell1};

const string kDSUName0 = "DSU0";
const string kDSUName1 = "DSU1";
const string kDSUName2 = "DSU2";
const string kDSUName3 = "DSU3";

const string kFile0 = "file0";
const string kFile1 = "file1";
const string kFile2 = "file2";
const string kFile3 = "file3";

void load_file(kmrnext::DataStore *ds, const string& file);
void print_ds(kmrnext::DataStore *ds, const kmrnext::View& v,
	      const string& name, int count);


//////////////////////////////////////////////////////////////////////////////
// Main starts from here.
//////////////////////////////////////////////////////////////////////////////
int
main(int argc, char **argv)
{
  kmrnext::KMRNext *next = kmrnext::KMRNext::init(argc, argv);

  /////////// Create DataStores
  kmrnext::DataStore *dsu0 = next->create_ds(kDimCell);
  kmrnext::DataStore *dsu1 = next->create_ds(kDimCell);
  kmrnext::DataStore *dsu2 = next->create_ds(kDimCell);
  kmrnext::DataStore *dsu3 = next->create_ds(kDimCell);
  dsu0->set(kDSCellSizes);
  dsu1->set(kDSCellSizes);
  dsu2->set(kDSCellSizes);
  dsu3->set(kDSCellSizes);
  if (kPrint) {
    cout << "0. Create DataStores" << endl;
    cout << "  DSCell0: " << dsu0->to_string() << endl;
    cout << "  DSCell1: " << dsu1->to_string() << endl;
    cout << "  DSCell2: " << dsu2->to_string() << endl;
    cout << "  DSCell3: " << dsu3->to_string() << endl;
    cout << endl;
  }

  /////////// Load data from files
  load_file(dsu0, kFile0);
  load_file(dsu1, kFile1);
  load_file(dsu2, kFile2);
  load_file(dsu3, kFile3);
  if (kPrint) {
    cout << "1. Load data to DataStores" << endl;
    kmrnext::View v(kDimCell);
    bool flag[kDimCell] = {false, false};
    v.set(flag);
    print_ds(dsu0, v, kDSUName0, 10);
    print_ds(dsu1, v, kDSUName1, 10);
    print_ds(dsu2, v, kDSUName2, 10);
    print_ds(dsu3, v, kDSUName3, 10);
    cout << endl;
  }

  /////////// Conmbine DataStores to create a high-dimensinal DataStore
  // 1. Test one DataStore
  kmrnext::DataStore *ds0 = next->create_ds(kDimCell + 1);
  vector<kmrnext::DataStore*> dsus;
  dsus.push_back(dsu0);
  ds0->set_from(dsus);
  if (kPrint) {
    cout << "2. Increase dimension by set_from()" << endl;
    kmrnext::View v(kDimCell + 1);
    bool flag[kDimCell + 1] = {false, false, false};
    v.set(flag);
    print_ds(ds0, v, "DS_Single", -1);
    cout << endl;
  }
  delete ds0;
  // 2. Test four DataStores
  kmrnext::DataStore *ds1 = next->create_ds(kDimCell + 1);
  dsus.push_back(dsu1);
  dsus.push_back(dsu2);
  dsus.push_back(dsu3);
  ds1->set_from(dsus);
  if (kPrint) {
    cout << "3. Combine four DataStores by set_from()" << endl;
    kmrnext::View v(kDimCell + 1);
    bool flag[kDimCell + 1] = {false, false, false};
    v.set(flag);
    print_ds(ds1, v, "DS_Four", -1);
    cout << endl;
  }
  delete ds1;

  /////////// Split a DataStore to low-dimensional DataStores
  kmrnext::DataStore *dsm0 = next->create_ds(kDimCell - 1);
  kmrnext::DataStore *dsm1 = next->create_ds(kDimCell - 1);
  kmrnext::DataStore *dsm2 = next->create_ds(kDimCell - 1);
  kmrnext::DataStore *dsm3 = next->create_ds(kDimCell - 1);
  kmrnext::DataStore *dsm4 = next->create_ds(kDimCell - 1);
  kmrnext::DataStore *dsm5 = next->create_ds(kDimCell - 1);
  kmrnext::DataStore *dsm6 = next->create_ds(kDimCell - 1);
  kmrnext::DataStore *dsm7 = next->create_ds(kDimCell - 1);
  kmrnext::DataStore *dsm8 = next->create_ds(kDimCell - 1);
  kmrnext::DataStore *dsm9 = next->create_ds(kDimCell - 1);
  vector<kmrnext::DataStore*> dsms;
  dsms.push_back(dsm0);
  dsms.push_back(dsm1);
  dsms.push_back(dsm2);
  dsms.push_back(dsm3);
  dsms.push_back(dsm4);
  dsms.push_back(dsm5);
  dsms.push_back(dsm6);
  dsms.push_back(dsm7);
  dsms.push_back(dsm8);
  dsms.push_back(dsm9);
  dsu0->split_to(dsms);
  if (kPrint) {
    cout << "4. Split a DataStore to 10 DataStores by split_to()"  << endl;
    kmrnext::View v(kDimCell - 1);
    bool flag[kDimCell - 1] = {false};
    v.set(flag);
    print_ds(dsms.at(0), v, "DS_Smallest", -1);
    print_ds(dsms.at(1), v, "DS_Smallest", -1);
    cout << endl;
  }
  delete dsm0;
  delete dsm1;
  delete dsm2;
  delete dsm3;
  delete dsm4;
  delete dsm5;
  delete dsm6;
  delete dsm7;
  delete dsm8;
  delete dsm9;

  delete dsu0;
  delete dsu1;
  delete dsu2;
  delete dsu3;
  kmrnext::KMRNext::finalize();
  return 0;
}


class DataLoader : public kmrnext::DataStore::Loader<string> {
public:
  int operator()(kmrnext::DataStore *ds, const string& file)
  {
    kmrnext::Key key(kDimCell);
    for (size_t i = 0; i < kDimCell0; i++) {
      key.set_dim(0, i);
      for (size_t j = 0; j < kDimCell1; j++) {
    	key.set_dim(1, j);
	long val;
	if (file == kFile0) {
	  val = 1;
	} else if (file == kFile1) {
	  val = 2;
	} else if (file == kFile2) {
	  val = 3;
	} else if (file == kFile3) {
	  val = 4;
	} else {
	  val = 0;
	}
	kmrnext::Data d(&val, sizeof(long));
	ds->add(key, d);
      }
    }
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
		 kmrnext::Key& key, vector<kmrnext::DataPack>& dps)
  {
    //cout << _padding << "Key: " << key.to_string() << endl;
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

void load_file(kmrnext::DataStore *ds, const string& file)
{
  vector<string> files;
  files.push_back(file);
  DataLoader loader;
  ds->load_files(files, loader);
}

void print_ds(kmrnext::DataStore *ds, const kmrnext::View& v,
	      const string& name, int count)
{
  cout << "  " << name << endl;
  DataStorePrinter dsp(count, "    ");
  ds->map(NULL, dsp, v);
}
