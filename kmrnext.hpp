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

    size_t size()
    {
      return _size;
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

    void copy_deep(const Data& src);

    // It returns a pointer to the stored data.
    void *value() { return _value; }

    // It returns size of the stored data.
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

    // It sets size of each dimension.
    virtual void set(size_t *val);

    // It adds a data to this DataStore.
    void add(const Key& key, const Data& data);

    // It gets a specified data from this DataStore.
    DataPack get(const Key& key);

    // It gets data whose keys are same when the specified view is applied.
    vector<DataPack>* get(const View& view, const Key& key);

    template <typename Mapper>
    void map(DataStore* outds, Mapper m, const View& view)
    {
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

      for (int i = 0; i < dpgroups.size(); i++) {
	vector<DataPack> &dps = dpgroups.at(i);
	Key viewed_key = key_to_viewed_key(dps.at(0).key, view);
	m(this, outds, viewed_key, dps);
      }
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
