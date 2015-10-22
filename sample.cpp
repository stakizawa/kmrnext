#include <iostream>
#include <sstream>
#include <vector>
//#include <typeinfo>

using namespace std;

template <std::size_t Dim>
class Key {
public:
  size_t value[Dim];

  string to_string() {
    ostringstream os;
    os << '<';
    for (std::size_t i = 0; i < Dimension; i++) {
      os << value[i];
      if (i < Dimension - 1) {
	os << ',';
      }
    }
    os << '>';
    return os.str();
  }

  static const std::size_t Dimension = Dim;
};

class Data {
public:
  void        *value;
  std::size_t value_size;
};

template <std::size_t Dim>
class View {
  bool _value[Dim];

public:
  View(bool *flags)
  {
    for (std::size_t i = 0; i < Dimension; i++) {
      _value[i] = flags[i];
    }
  }

  string to_string() {
    ostringstream os;
    os << '<';
    for (std::size_t i = 0; i < Dimension; i++) {
      os << _value[i];
      if (i < Dimension - 1) {
	os << ',';
      }
    }
    os << '>';
    return os.str();
  }

  static const std::size_t Dimension = Dim;
};

template <std::size_t Dim>
class DataStore {
public:
  DataStore(std::size_t *sizes)
    : _data(NULL), _data_size(1) {
    for (std::size_t i = 0; i < Dimension; i++) {
      _dim_sizes[i] = sizes[i];
      _data_size *= sizes[i];
    }
    _data = static_cast<Data*>(malloc(sizeof(Data) * _data_size));
  }

  void add(Key<Dim>& key, void *data, std::size_t data_size)
  {
    size_t idx = calc_index(key);
    Data *d = &(_data[idx]);
    d->value = static_cast<void*>(malloc(data_size));
    memcpy(d->value, data, data_size);
    d->value_size = data_size;
  }

  Data& get(Key<Dim>& key)
  {
    size_t idx = calc_index(key);
    return _data[idx];
  }

  template <class Functor>
  void load_files(const std::vector<string>& files, Functor f)
  {
    for (size_t i = 0; i < files.size(); i++) {
      string file = files[i];
      f(this, file);
    }
  }

  string to_string()
  {
    ostringstream os;
    os << '<';
    for (std::size_t i = 0; i < Dimension; i++) {
      os << _dim_sizes[i];
      if (i < Dimension - 1) {
	os << ',';
      }
    }
    os << '>';
    return os.str();
  }

  static const std::size_t Dimension = Dim;

private:
  std::size_t _dim_sizes[Dim];
  Data *_data;
  size_t _data_size;

  size_t calc_index(Key<Dim>& key)
  {
    size_t idx = 0;
    for (size_t i = 0; i < Dimension; i++) {
      size_t offset = 1;
      for (size_t j = i+1; j < Dimension; j++) {
	offset *= _dim_sizes[j];
      }
      idx += key.value[i] * offset;
    }
    return idx;
  }
};


///////////////////////////////////////////////////////////////////////////////
// Main starts from here.
///////////////////////////////////////////////////////////////////////////////
const int Dimension = 3;
typedef Key<Dimension> Key3;
typedef DataStore<Dimension> DS3;
typedef View<Dimension> V3;

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
	  ds->add(key, &val, sizeof(long));
	}
      }
    }
    return 0;
  }
};

int
main()
{
  cout << DS3::Dimension << endl;

  ///////////  Create a DataStore
  size_t sizes[Dimension] = {Dim0, Dim1, Dim2};
  DS3 ds1(sizes);
  cout << ds1.to_string() << endl;

  ///////////  Load data contents from a file
  std::vector<string> files;
  files.push_back("dummy1");
  //  files.push_back("dummy2");
  Loader loader;
  ds1.load_files(files, loader);

  ///////////  Get a data from a DataStore
  Key3 key;
  key.value[0] = 2; key.value[1] = 2; key.value[2] = 2;
  Data d1 = ds1.get(key);
  cout << "Value: " << *(long *)d1.value << endl;
  key.value[0] = 2; key.value[1] = 2; key.value[2] = 3;
  d1 = ds1.get(key);
  cout << "Value: " << *(long *)d1.value << endl;

  ///////////  Set a View then map function
  bool flags1[Dimension] = {true, true, true};
  V3 v1(flags1);
  cout << v1.to_string() << endl;
  //  ds1.set_view(v1);
  //  ds1.map(mapper1);

  bool flags2[Dimension] = {true, false, false};
  V3 v2(flags2);
  // ds1.set_view(v2);
  // ds1.map(mapper2);

  return 0;
}
