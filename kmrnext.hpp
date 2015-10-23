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
    Dimensional(size_t size) : _size(size) {}

    virtual void set(T *val)
    {
      for (size_t i = 0; i < _size; i++) {
	_value[i] = val[i];
      }
    }

    T dim(size_t idx) const
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

  ///////////////////////////////////////////////////////////////////////////
  // Key class
  ///////////////////////////////////////////////////////////////////////////
  class Key : public Dimensional<size_t> {
  public:
    Key(size_t size) : Dimensional<size_t>(size) {}
  };

  ///////////////////////////////////////////////////////////////////////////
  // A class that define view of data
  ///////////////////////////////////////////////////////////////////////////
  class View : public Dimensional<bool> {
  public:
    View(size_t size) : Dimensional<bool>(size) {}
  };

  ///////////////////////////////////////////////////////////////////////////
  // A class that define view of data
  ///////////////////////////////////////////////////////////////////////////
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

  ///////////////////////////////////////////////////////////////////////////
  // A class that stores a data ans its key
  ///////////////////////////////////////////////////////////////////////////
  class DataPack {
  public:
    Key key;
    Data *data;

    DataPack(const Key k, Data *d) : key(k), data(d) {}
  };

  class DataStore : public Dimensional<size_t> {
  public:
    DataStore(size_t size)
      : Dimensional<size_t>(size), _data(NULL), _data_size(0) {}

    virtual void set(size_t *val)
    {
      _data_size = 1;
      for (size_t i = 0; i < _size; i++) {
	_value[i] = val[i];
	_data_size *= val[i];
      }
      _data = static_cast<Data*>(malloc(sizeof(Data) * _data_size));
    }

    void add(const Key& key, const Data& data)
    {
      size_t idx = key_to_index(key);
      Data *d = &(_data[idx]);
      d->copy_deep(data);
    }

    DataPack get(const Key& key)
    {
      size_t idx = key_to_index(key);
      return DataPack(key, &(_data[idx]));
    }

    vector<DataPack>* get(const View& view, const Key& key)
    {
      vector<DataPack> *dps = new vector<DataPack>();
      for (size_t i = 0; i < _data_size; i++) {
	Key tmpkey = index_to_key(i);
	bool store = true;
	for (size_t j = 0; j < _size; j++) {
	  if (view.dim(j) && key.dim(j) != tmpkey.dim(j)) {
	    store = false;
	    break;
	  }
	}
	if (store) {
	  //	  dps->push_back(_data[i]);
	  dps->push_back(DataPack(tmpkey, &(_data[i])));
	}
      }
      return dps;
    }

    template <typename Mapper>
    void map(DataStore& outds, Mapper m, const View& view)
    {
    }

    template <typename Loader>
    void load_files(const vector<string>& files, Loader f)
    {
      for (size_t i = 0; i < files.size(); i++) {
	string file = files[i];
	f(this, file);
      }
    }

  private:
    Data *_data;
    size_t _data_size;

    size_t key_to_index(const Key& key)
    {
      size_t idx = 0;
      for (size_t i = 0; i < _size; i++) {
	size_t offset = 1;
	for (size_t j = i+1; j < _size; j++) {
	  offset *= _value[j];
	}
	idx += key.dim(i) * offset;
      }
      return idx;
    }

    Key index_to_key(const size_t index)
    {
      Key key(_size);
      size_t _index = index;
      for (size_t i = 0; i < _size; i++) {
	size_t length = 1;
	for (size_t j = i+1; j < _size; j++) {
	  length *= _value[j];
	}
	key.set_dim(i, _index / length);
	_index %= length;
      }
      return key;
    }
  };

}

#endif
