#include <iostream>
#include <sstream>
#include <vector>

using namespace std;

template <size_t Dim>
class Key {
public:
  size_t value[Dim];

  Key() {}

  Key(size_t *val)
  {
    for (size_t i = 0; i < Dimension; i++) {
      value[i] = val[i];
    }
  }

  string to_string()
  {
    ostringstream os;
    os << '<';
    for (size_t i = 0; i < Dimension; i++) {
      os << value[i];
      if (i < Dimension - 1) {
	os << ',';
      }
    }
    os << '>';
    return os.str();
  }

  static const size_t Dimension = Dim;
  static const size_t Any = -1;
};

class Data {
  void *_value;
  size_t _value_size;

public:
  Data(void *value, const size_t value_size)
    : _value(value), _value_size(value_size) {}

  void copy_deep(const Data& src)
  {
    _value = static_cast<void*>(malloc(src._value_size));
    memcpy(_value, src._value, src._value_size);
    _value_size = src._value_size;
  }

  void *value() { return _value; }

  size_t size() { return _value_size; }
};

template <size_t Dim>
class View {
  bool _value[Dim];

public:
  View(bool *flags)
  {
    for (size_t i = 0; i < Dimension; i++) {
      _value[i] = flags[i];
    }
  }

  bool dimension(int dim) const
  {
    return _value[dim];
  }

  string to_string()
  {
    ostringstream os;
    os << '<';
    for (size_t i = 0; i < Dimension; i++) {
      os << _value[i];
      if (i < Dimension - 1) {
	os << ',';
      }
    }
    os << '>';
    return os.str();
  }

  static const size_t Dimension = Dim;
};

template <size_t Dim>
class DataStore {
public:
  DataStore(size_t *sizes)
    : _data(NULL), _data_size(1)
  {
    for (size_t i = 0; i < Dimension; i++) {
      _dim_sizes[i] = sizes[i];
      _data_size *= sizes[i];
    }
    _data = static_cast<Data*>(malloc(sizeof(Data) * _data_size));
  }

  void add(const Key<Dim>& key, const Data& data)
  {
    size_t idx = key_to_index(key);
    Data *d = &(_data[idx]);
    d->copy_deep(data);
  }

  Data& get(const Key<Dim>& key)
  {
    size_t idx = key_to_index(key);
    return _data[idx];
  }

  vector<Data>* get(const View<Dim>& view, const Key<Dim>& key)
  {
    vector<Data> *dvec = new vector<Data>();
    for (size_t i = 0; i < _data_size; i++) {
      Key<Dim> *k = index_to_key(i);
      bool store = true;
      for (size_t j = 0; j < Dimension; j++) {
	if (view.dimension(j) && key.value[j] != k->value[j]) {
	  store = false;
	  break;
	}
      }
      delete k;
      if (store) {
	dvec->push_back(_data[i]);
      }
    }
    return dvec;
  }

  template <typename Loader>
  void load_files(const vector<string>& files, Loader f)
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
    for (size_t i = 0; i < Dimension; i++) {
      os << _dim_sizes[i];
      if (i < Dimension - 1) {
	os << ',';
      }
    }
    os << '>';
    return os.str();
  }

  static const size_t Dimension = Dim;

private:
  size_t _dim_sizes[Dim];
  Data *_data;
  size_t _data_size;

  size_t key_to_index(const Key<Dim>& key)
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

  Key<Dim>* index_to_key(const size_t index)
  {
    Key<Dim> *key = new Key<Dim>();
    size_t _index = index;
    for (size_t i = 0; i < Dimension; i++) {
      size_t length = 1;
      for (size_t j = i+1; j < Dimension; j++) {
	length *= _dim_sizes[j];
      }
      key->value[i] = _index / length;
      _index %= length;
    }
    return key;
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
	  Data d(&val, sizeof(long));
	  ds->add(key, d);
	}
      }
    }
    return 0;
  }
};

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
  Data d1 = ds1.get(key1);
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
  vector<Data> *dvec1 = ds1.get(v1, key1);
  cout << "dvec1" << endl;
  cout << "  size: " << dvec1->size() << endl;
  // TODO print

  ///////////  Set a View then map function
  cout << v1.to_string() << endl;
  //  ds1.set_view(v1);
  //  ds1.map(v1, mapper1);

  // ds1.set_view(v2);
  // ds1.map(mapper2);

  return 0;
}
