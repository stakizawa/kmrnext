/// \file
/// Example of matrix transpose
///
/// It transposes a matrix.
/// In case of SERIAL backend, the matrix size is constant and 4x4.
/// In case of KMR backend, the matrix size is nprocs x nprocs, where
/// nprocs is the number of nodes.
#include <iostream>
#include <cassert>
#include "kmrnext.hpp"

using namespace std;
using namespace kmrnext;

const int kDimMatrix = 2;
const size_t kMatrixSizeInSerial = 4;

int rank = 0;
size_t matrix_size = 0;

void load_data(DataStore* ds);
void transpose(DataStore* ids, DataStore* ods);
void print_matrix(DataStore* ds);
void print_line(string& str);
void print_line(const char *str);

//////////////////////////////////////////////////////////////////////////////
// Main starts from here.
//////////////////////////////////////////////////////////////////////////////
int
main(int argc, char **argv)
{
  KMRNext *next = KMRNext::init(argc, argv);
#ifdef BACKEND_SERIAL
  matrix_size = kMatrixSizeInSerial;
#endif
#ifdef BACKEND_KMR
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, (int*)&matrix_size);
#endif

  DataStore* ds0 = next->create_ds(kDimMatrix);
  size_t msize[kDimMatrix] = {matrix_size, matrix_size};
  ds0->set(msize);
  DataStore* ds1 = next->create_ds(kDimMatrix);
  ds1->set(msize);
  DataStore* ds2 = next->create_ds(kDimMatrix);
  ds2->set(msize);

  // load data to ds0
  load_data(ds0);
  print_line("The original data.");
  print_matrix(ds0);

  // transpose ds0 and write outputs to ds1
  transpose(ds0, ds1);
  print_line("Transpose once.");
  print_matrix(ds1);
  delete ds0;

  // transpose ds1 and write outputs to ds2
  transpose(ds1, ds2);
  print_line("Transpose twice.");
  print_matrix(ds2);
  delete ds1;
  delete ds2;

  KMRNext::finalize();
  return 0;
}


class DataLoader : public DataStore::Loader<long> {
public:
  int operator()(DataStore* ds, const long& num)
  {
    int val = (int)num + 1;
    Data data((void*)&val, sizeof(int));

    Key key(kDimMatrix);
    key.set_dim(0, num);
    for (size_t i = 0; i < matrix_size; i++) {
      key.set_dim(1, i);
      ds->add(key, data);
    }
    return 0;
  }
};

void load_data(DataStore* ds)
{
  vector<long> dlist;
  for (size_t i = 0; i < matrix_size; i++) {
    dlist.push_back((int)i);
  }
  DataLoader dl;
  ds->load_integers(dlist, dl);
}

class Transposer : public DataStore::Mapper {
public:
  int operator()(DataStore* inds, DataStore* outds,
		 Key& key, vector<DataPack>& dps,
		 DataStore::MapEnvironment& env)
  {
#if DEBUG
    assert(dps.size() == 1);
#endif

    for (vector<DataPack>::iterator itr = dps.begin(); itr != dps.end();
	 itr++) {
      Key ikey = itr->key();
      Key okey(kDimMatrix);
      okey.set_dim(0, ikey.dim(1));
      okey.set_dim(1, ikey.dim(0));
      Data data(itr->data()->value(), itr->data()->size());
      outds->add(okey, data);
    }
    return 0;
  }
};

void transpose(DataStore* ids, DataStore* ods)
{
  Transposer mapper;
  View view(kDimMatrix);
  bool view_flag[kDimMatrix] = {true, true};
  view.set(view_flag);
  ids->map(mapper, view, ods);
}

void print_matrix(DataStore* ds)
{
#if DEBUG
  ostringstream os;
  Key key(kDimMatrix);
  for (size_t i = 0; i < matrix_size; i++) {
    key.set_dim(0, i);
    for (size_t j = 0; j < matrix_size; j++) {
      key.set_dim(1, j);
      DataPack dp = ds->get(key);
      Data *data = dp.data();
#ifdef BACKEND_SERIAL
      os << *(int*)data->value() << " ";
#elif defined BACKEND_KMR
      os << *(int*)data->value() << "(" << data->owner() << ") ";
#endif
    }
    os << endl;
  }

#ifdef BACKEND_KMR
  if (rank == 0) {
    cout << "Values inside parentheses is the rank of data owner." << endl;
  }
#endif
  string str0 = os.str();
  print_line(str0);
#endif
  return;
}

void print_line(string& str)
{
  if (rank == 0) {
    cout << str << endl;
  }
}

void print_line(const char *str)
{
  string str0(str);
  print_line(str0);
}
