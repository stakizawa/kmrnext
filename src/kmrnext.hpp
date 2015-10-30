#ifndef KMRNEXT_HPP
#define KMRNEXT_HPP

#include <iostream>
#include <sstream>
#include <stdexcept>
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
    explicit Dimensional(size_t siz) : _size(siz)
    {
      if (_size > MaxDimensionSize) {
	ostringstream os;
	os << "Dimension size should be less than " << (MaxDimensionSize + 1);
	throw runtime_error(os.str());
      }
    }

    virtual void set(const T *val)
    {
      for (size_t i = 0; i < _size; i++) {
	_value[i] = val[i];
      }
    }

    size_t size()
    {
      return _size;
    }

    T dim(size_t idx) const
    {
      if (idx >= _size) {
	ostringstream os;
	os << "Index should be less than dimension size: " << _size;
	throw runtime_error(os.str());
      }
      return _value[idx];
    }

    void set_dim(size_t idx, T val)
    {
      if (idx >= _size) {
	ostringstream os;
	os << "Index should be less than dimension size: " << _size;
	throw runtime_error(os.str());
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
    explicit Key(size_t siz) : Dimensional<size_t>(siz) {}
  };

  ///////////////////////////////////////////////////////////////////////////
  // A class that define view of data
  ///////////////////////////////////////////////////////////////////////////
  class View : public Dimensional<bool> {
  public:
    explicit View(size_t siz) : Dimensional<bool>(siz) {}
  };

  ///////////////////////////////////////////////////////////////////////////
  // A class that define view of data
  ///////////////////////////////////////////////////////////////////////////
  class Data {
    void *_value;
    size_t _value_size;

  public:
    Data(void *val, const size_t val_siz)
      : _value(val), _value_size(val_siz), _value_allocated(false) {}

    ~Data();

    void copy_deep(const Data& src);

    // It returns a pointer to the stored data.
    void *value() { return _value; }

    // It returns size of the stored data.
    size_t size() { return _value_size; }

  private:
    bool _value_allocated;
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
    explicit DataStore(size_t siz)
      : Dimensional<size_t>(siz), _data(NULL), _data_size(0),
      _data_allocated(false) {}

    ~DataStore();

    // It sets size of each dimension.
    virtual void set(const size_t *val);

    // It adds a data to this DataStore.
    void add(const Key& key, const Data& data);

    // It gets a specified data from this DataStore.
    DataPack get(const Key& key);

    // It gets data whose keys are same when the specified view is applied.
    vector<DataPack>* get(const View& view, const Key& key);

    // It sets Data from DataStores.
    void set_from(const vector<DataStore*>& dslist);

    // It splits DataStore to low-dimensional DataStores.
    void split_to(vector<DataStore*>& dslist);

    template <typename Mapper>
    void map(DataStore* outds, Mapper m, const View& view)
    {
      if (_data_size == 0) {
	return;
      }

      size_t nkeys = 1;
      for (size_t i = 0; i < _size; i++) {
	if (view.dim(i)) {
	  nkeys *= _value[i];
	}
      }

      vector< vector<DataPack> > dpgroups(nkeys);

      for (size_t i = 0; i < _data_size; i++) {
	Key tmpkey = index_to_key(i);
	size_t viewed_idx = key_to_viwed_index(tmpkey, view);
	vector<DataPack>& dps = dpgroups.at(viewed_idx);
	dps.push_back(DataPack(tmpkey, &(_data[i])));
      }

      for (size_t i = 0; i < dpgroups.size(); i++) {
	vector<DataPack> &dps = dpgroups.at(i);
	Key viewed_key = key_to_viewed_key(dps.at(0).key, view);
	m(this, outds, viewed_key, dps);
      }
    }

    template <typename Loader>
    void load_files(const vector<string>& files, Loader f)
    {
      load_array(files, f);
      // for (size_t i = 0; i < files.size(); i++) {
      // 	f(this, files[i]);
      // }
    }

    template <typename Type, typename Loader>
    void load_array(const vector<Type>& array, Loader f)
    {
      if (array.size() == 1) {
	f(this, array[0]);
	return;
      }

      // Check if the size of array is same as the multiple of dimension.
      size_t mi;
      size_t mdim = 1;
      for (mi = 0; mi < _size; mi++) {
	mdim *= _value[mi];
	if (mdim == array.size()) {
	  break;
	}
      }
      if (mdim != array.size()) {
	throw runtime_error("The size of array should be match the "
			    "size of dimension.");
      }

      // Calculate size and dimension of Sub DS
      size_t sub_ds_dim = _size - (mi + 1);
      size_t sub_ds_dims[MaxDimensionSize];
      size_t sub_ds_siz = 0;
      for (size_t i = 0; i < sub_ds_dim; i++) {
	sub_ds_dims[i] = _value[i + mi + 1];
	sub_ds_siz += _value[i + mi + 1];
      }

      size_t offset = 0;
      for (size_t i = 0; i < array.size(); i++) {
      	DataStore ds(sub_ds_dim);
      	ds.set(sub_ds_dims, this->_data + offset);
      	f(&ds, array[i]);
	offset += sub_ds_siz;
      }
    }

  private:
    Data *_data;
    size_t _data_size;
    bool _data_allocated;

    // It sets size of each dimension.
    // It just sets pointer to data, not performs malloc and memcpy.
    void set(const size_t *val, Data *dat_ptr);

    // It returns the index of Data calculated from the specified Key.
    size_t key_to_index(const Key& key);

    // It returns a Key of the specified indexed Data.
    Key index_to_key(const size_t index);

    // It returns the index of Data calculated from the specified Key when
    // the specified view is applied.
    size_t key_to_viwed_index(const Key& key, const View& view);

    // It converts the specified key by applying the specified View.
    Key key_to_viewed_key(const Key& key, const View& view);
  };

}

#endif