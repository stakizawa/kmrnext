#include <sstream>
#include <stdexcept>
#include <cstdlib>
#include <cstring>
#include "kmrnext.hpp"

namespace Next {
  Data::~Data()
  {
    if (_value_allocated) {
      free(_value);
    }
  }

  void Data::copy_deep(const Data& src)
  {
    _value = static_cast<void*>(malloc(src._value_size));
    memcpy(_value, src._value, src._value_size);
    _value_size = src._value_size;
  }

  DataStore::~DataStore()
  {
    if (_data_allocated) {
      free(_data);
    }
  }

  void DataStore::set(const size_t *val)
  {
    if (_data_size != 0) {
      throw runtime_error("DataStore is already initialized.");
    }

    _data_size = 1;
    for (size_t i = 0; i < _size; i++) {
      _value[i] = val[i];
      _data_size *= val[i];
    }
    _data = static_cast<Data*>(malloc(sizeof(Data) * _data_size));
    _data_allocated = true;
  }

  void DataStore::add(const Key& key, const Data& data)
  {
    size_t idx = key_to_index(key);
    Data *d = &(_data[idx]);
    d->copy_deep(data);
  }

  DataPack DataStore::get(const Key& key)
  {
    size_t idx = key_to_index(key);
    return DataPack(key, &(_data[idx]));
  }

  vector<DataPack>* DataStore::get(const View& view, const Key& key)
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
	dps->push_back(DataPack(tmpkey, &(_data[i])));
      }
    }
    return dps;
  }

  void DataStore::set_from(const vector<DataStore*>& dslist)
  {
    if (dslist.size() == 0) {
      throw runtime_error("There should be at least one DataStore.");
    }
    if (_data_size != 0) {
      throw runtime_error("DataStore is already initialized.");
    }
    // TODO check if the size of each DataStore in dslist is same.

    _value[0] = dslist.size();
    DataStore *ds0 = dslist.at(0);
    for (size_t i = 1; i < _size; i++) {
      _value[i] = ds0->_value[i-1];
    }

    _data_size = 1;
    for (size_t i = 0; i < _size; i++) {
      _data_size *= _value[i];
    }
    _data = static_cast<Data*>(malloc(sizeof(Data) * _data_size));

    size_t offset = 0;
    for (size_t i = 0; i < dslist.size(); i++) {
      DataStore *src = dslist.at(i);
      memcpy(_data + offset, src->_data, sizeof(Data) * src->_data_size);
      offset += src->_data_size;
    }
  }

  void DataStore::split_to(vector<DataStore*>& dslist)
  {
    if (_data_size == 0) {
      throw runtime_error("Data should be set.");
    }
    if (_size < 2) {
      throw runtime_error("DataStore can't be split.");
    }
    if (_value[0] != dslist.size()) {
      ostringstream os;
      os << "DataStore vector size should be " << _value[0]
	 << ", but " << dslist.size() << ".";
      throw runtime_error(os.str());
    }
    // TODO check if the size of each DataStore in dslist is same.

    size_t split_dims[MaxDimensionSize];
    for (size_t i = 1; i < _size; i++) {
      split_dims[i-1] = _value[i];
    }

    size_t offset = 0;
    for (size_t i = 0; i < dslist.size(); i++) {
      DataStore *dst = dslist.at(i);
      dst->set(split_dims);
      memcpy(dst->_data, _data + offset, sizeof(Data) * dst->_data_size);
      offset += dst->_data_size;
    }
  }

  void DataStore::set(const size_t *val, Data *dat_ptr)
  {
    if (_data_size != 0) {
      throw runtime_error("DataStore is already initialized.");
    }

    _data_size = 1;
    for (size_t i = 0; i < _size; i++) {
      _value[i] = val[i];
      _data_size *= val[i];
    }
    _data = dat_ptr;
  }

  size_t DataStore::key_to_index(const Key& key)
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

  Key DataStore::index_to_key(const size_t index)
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

  size_t DataStore::key_to_viwed_index(const Key& key, const View& view)
  {
    size_t idx = 0;
    for (size_t i = 0; i < _size; i++) {
      if (view.dim(i)) {
	size_t offset = 1;
	for (size_t j = i+1; j < _size; j++) {
	  if (view.dim(j)) {
	    offset *= _value[j];
	  }
	}
	idx += key.dim(i) * offset;
      }
    }
    return idx;
  }

  Key DataStore::key_to_viewed_key(const Key& key, const View& view)
  {
    size_t viewed_key_size = 0;
    for (size_t i = 0; i < _size; i++) {
      if (view.dim(i)) {
	viewed_key_size += 1;
      }
    }

    Key viewed_key(viewed_key_size);
    size_t viewed_key_idx = 0;
    for (size_t i = 0; i < _size; i++) {
      if (view.dim(i)) {
	viewed_key.set_dim(viewed_key_idx, key.dim(i));
	viewed_key_idx += 1;
      }
    }
    return viewed_key;
  }

#if 0
  // A version where dimension of the viewed key is same as the source key
  Key DataStore::key_to_viewed_key(const Key& key, const View& view)
  {
    Key viewed_key(_size);
    for (size_t i = 0; i < _size; i++) {
      if (view.dim(i)) {
	viewed_key.set_dim(i, key.dim(i));
      } else {
	viewed_key.set_dim(i, 0);
      }
    }
    return viewed_key;
  }
#endif

}
