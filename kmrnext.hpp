#ifndef KMRNEXT_HPP
#define KMRNEXT_HPP

#include <iostream>
#include <sstream>
#include <vector>

namespace Next {
  using namespace std;

  const size_t MaxDimensionSize = 8;

  template <typename T>
  class Dimensional {
  protected:
    size_t _size;
    T _value[MaxDimensionSize];

  public:
    T dim(size_t idx) const // TODO rename to operator[]?
    {
      if (idx >= _size) {
	// TODO through exception
      }
      return _value[idx];
    }

    void set_dim(size_t idx, T val)
    {
      if (idx >= _size) {
	// TODO through exception
      }
      _value[idx] = val;
    }

    string to_string()
    {
      ostringstream os;
      os << '<';
      for (size_t i = 0; i < _size; i++) {
	os << _value[i];
	if (i < _size - 1) {
	  os << ',';
	}
      }
      os << '>';
      return os.str();
    }
  };

  class Key : public Dimensional<size_t> {
  public:
    Key(size_t dim_size)
    {
      _size = dim_size;
    }

    void set(size_t *val)
    {
      for (size_t i = 0; i < _size; i++) {
	_value[i] = val[i];
      }
    }
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

    void add(const Key& key, const Data& data)
    {
      size_t idx = key_to_index(key);
      Data *d = &(_data[idx]);
      d->copy_deep(data);
    }

    Data& get(const Key& key)
    {
      size_t idx = key_to_index(key);
      return _data[idx];
    }

    vector<Data>* get(const View<Dim>& view, const Key& key)
    {
      vector<Data> *dvec = new vector<Data>();
      for (size_t i = 0; i < _data_size; i++) {
	Key *k = index_to_key(i);
	bool store = true;
	for (size_t j = 0; j < Dimension; j++) {
	  if (view.dimension(j) && key.dim(j) != k->dim(j)) {
	    store = false;
	    break;
	  }
	}
	if (store) {
	  dvec->push_back(_data[i]);
	}
	delete k;
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

    size_t key_to_index(const Key& key)
    {
      size_t idx = 0;
      for (size_t i = 0; i < Dimension; i++) {
	size_t offset = 1;
	for (size_t j = i+1; j < Dimension; j++) {
	  offset *= _dim_sizes[j];
	}
	idx += key.dim(i) * offset;
      }
      return idx;
    }

    Key* index_to_key(const size_t index)
    {
      Key *key = new Key(Dimension);
      size_t _index = index;
      for (size_t i = 0; i < Dimension; i++) {
	size_t length = 1;
	for (size_t j = i+1; j < Dimension; j++) {
	  length *= _dim_sizes[j];
	}
	key->set_dim(i, _index / length);
	_index %= length;
      }
      return key;
    }
  };
}

#endif
