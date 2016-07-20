#include <iostream>
#include "kmrnext.hpp"
#ifdef BACKEND_KMR
#include <mpi.h>
#endif

using namespace std;

int rank = 0;

const bool kPrint = true;
const int  kDumpCount = 20;

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

void load_file(kmrnext::DataStore* ds, const string& file);
void print_line(string str);
void print_data_store(kmrnext::DataStore* ds, string name, string padding,
		      int count);


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
    print_line("0. Create DataStores");
    print_line("  DSCell0: " + dsu0->to_string());
    print_line("  DSCell1: " + dsu1->to_string());
    print_line("  DSCell2: " + dsu2->to_string());
    print_line("  DSCell3: " + dsu3->to_string());
    print_line("");
  }

  /////////// Load data from files
  load_file(dsu0, kFile0);
  load_file(dsu1, kFile1);
  load_file(dsu2, kFile2);
  load_file(dsu3, kFile3);
  if (kPrint) {
    print_line("1. Load data to DataStores");
    print_data_store(dsu0, kDSUName0, "  ", kDumpCount);
    print_data_store(dsu1, kDSUName1, "  ", kDumpCount);
    print_data_store(dsu2, kDSUName2, "  ", kDumpCount);
    print_data_store(dsu3, kDSUName3, "  ", kDumpCount);
    print_line("");
  }

  /////////// Conmbine DataStores to create a high-dimensinal DataStore
  // 1. Test one DataStore
  kmrnext::DataStore *ds0 = next->create_ds(kDimCell + 1);
  vector<kmrnext::DataStore*> dsus;
  dsus.push_back(dsu0);
  ds0->set_from(dsus);
  if (kPrint) {
    print_line("2. Increase dimension by set_from()");
    print_data_store(ds0, "DS_Single", "  ", -1);
    print_line("");
  }
  delete ds0;
  // 2. Test four DataStores
  kmrnext::DataStore *ds1 = next->create_ds(kDimCell + 1);
  dsus.push_back(dsu1);
  dsus.push_back(dsu2);
  dsus.push_back(dsu3);
  ds1->set_from(dsus);
  if (kPrint) {
    print_line("3. Combine four DataStores by set_from()");
    print_data_store(ds1, "DS_Four", "  ", -1);
    print_line("");
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
    print_line("4. Split a DataStore to 10 DataStores by split_to()" );
    print_data_store(dsms.at(0), "DS_Smallest0", "  ", -1);
    print_data_store(dsms.at(1), "DS_Smallest1", "  ", -1);
    print_line("");
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
  int operator()(kmrnext::DataStore* ds, const string& file)
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

void load_file(kmrnext::DataStore* ds, const string& file)
{
  vector<string> files;
  files.push_back(file);
  DataLoader loader;
  ds->load_files(files, loader);
}

void print_line(string str) {
  if (rank != 0) return;
  cout << str << endl;
}

// A class for dumping a DataPack
class DPPrinter : public kmrnext::DataPack::Dumper {
  string padding_;
public:
  DPPrinter(string padding) : padding_(padding) {}
  string operator()(kmrnext::DataPack& dp) {
    ostringstream os;
    os << padding_ << dp.key().to_string() << " : "
       << *static_cast<long*>(dp.data()->value()) << endl;
    return os.str();
  }
};

void print_data_store(kmrnext::DataStore* ds, string name, string padding,
		      int count) {
  long ds_count = ds->count();
  DPPrinter printer(padding + padding);
  string ds_str = ds->dump(printer);
  if (rank != 0) return;

  cout << padding << name << ": " << ds->to_string() << endl;
  cout << padding << "  Count of data in the DataStore: " << ds_count << endl;
  if (count > 0) {
    cout << padding << "  Values (top" << count << ")" << endl;
  } else {
    cout << padding << "  Values (all)" << endl;
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
