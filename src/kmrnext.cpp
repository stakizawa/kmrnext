#include <sstream>
#include <stdexcept>
#include <cstdlib>
#include <cstring>
#include "kmrnext.hpp"

namespace kmrnext {
  Data::~Data()
  {
    if (value_allocated_) {
      free(value_);
    }
  }

  void Data::copy_deep(const Data& src)
  {
    if (value_ != NULL) {
      throw runtime_error("Data is already set value.");
    }
    value_ = static_cast<void*>(calloc(src.value_size_, sizeof(void*)));
    memcpy(value_, src.value_, src.value_size_);
    value_size_ = src.value_size_;
    value_allocated_ = true;
  }

  DataStore::~DataStore()
  {
    if (data_allocated_) {
      free(data_);
    }
  }

  void DataStore::set(const size_t *val)
  {
    if (data_size_ != 0) {
      throw runtime_error("DataStore is already initialized.");
    }

    data_size_ = 1;
    for (size_t i = 0; i < size_; i++) {
      value_[i] = val[i];
      data_size_ *= val[i];
    }
    data_ = static_cast<Data*>(calloc(data_size_, sizeof(Data)));
    data_allocated_ = true;
  }

  void DataStore::add(const Key& key, const Data& data)
  {
    check_key_range(key);
    size_t idx = key_to_index(key);
    Data *d = &(data_[idx]);
    d->copy_deep(data);
  }

  DataPack DataStore::get(const Key& key)
  {
    check_key_range(key);
    size_t idx = key_to_index(key);
    return DataPack(key, &(data_[idx]));
  }

  vector<DataPack>* DataStore::get(const View& view, const Key& key)
  {
    check_key_range(key);
    vector<DataPack> *dps = new vector<DataPack>();
    for (size_t i = 0; i < data_size_; i++) {
      Key tmpkey = index_to_key(i);
      bool store = true;
      for (size_t j = 0; j < size_; j++) {
	if (view.dim(j) && key.dim(j) != tmpkey.dim(j)) {
	  store = false;
	  break;
	}
      }
      if (store) {
	dps->push_back(DataPack(tmpkey, &(data_[i])));
      }
    }
    return dps;
  }

  void DataStore::set_from(const vector<DataStore*>& dslist)
  {
    if (dslist.size() == 0) {
      throw runtime_error("There should be at least one DataStore.");
    }
    if (data_size_ != 0) {
      throw runtime_error("DataStore is already initialized.");
    }
    {
      // Check each DataStore in the vector
      size_t expected_dim_size = size_ - 1;
      size_t expected_data_size = 0;
      for (size_t i = 0; i < dslist.size(); i++) {
	DataStore *src = dslist.at(i);
	if (expected_dim_size != src->size_) {
	  throw runtime_error("Dimension size of one of DataStore is wrong.");
	}
	size_t calc_data_size = 1;
	for (size_t j = 0; j < expected_dim_size; j++) {
	  calc_data_size *= src->value_[j];
	}
	if (i == 0) {
	  expected_data_size = calc_data_size;
	} else {
	  if (expected_data_size != calc_data_size) {
	    throw runtime_error("Data count one of DataStore is wrong.");
	  }
	}
      }
    }

    value_[0] = dslist.size();
    DataStore *ds0 = dslist.at(0);
    for (size_t i = 1; i < size_; i++) {
      value_[i] = ds0->value_[i-1];
    }

    data_size_ = 1;
    for (size_t i = 0; i < size_; i++) {
      data_size_ *= value_[i];
    }
    data_ = static_cast<Data*>(calloc(data_size_, sizeof(Data)));

    size_t offset = 0;
    for (size_t i = 0; i < dslist.size(); i++) {
      DataStore *src = dslist.at(i);
      memcpy(data_ + offset, src->data_, sizeof(Data) * src->data_size_);
      offset += src->data_size_;
    }
  }

  void DataStore::split_to(vector<DataStore*>& dslist)
  {
    if (data_size_ == 0) {
      throw runtime_error("Data should be set.");
    }
    if (size_ < 2) {
      throw runtime_error("DataStore can't be split.");
    }
    if (value_[0] != dslist.size()) {
      ostringstream os;
      os << "DataStore vector size should be " << value_[0]
	 << ", but " << dslist.size() << ".";
      throw runtime_error(os.str());
    }
    // TODO check if the size of each DataStore in dslist is same.

    size_t split_dims[kMaxDimensionSize];
    for (size_t i = 1; i < size_; i++) {
      split_dims[i-1] = value_[i];
    }

    size_t offset = 0;
    for (size_t i = 0; i < dslist.size(); i++) {
      DataStore *dst = dslist.at(i);
      dst->set(split_dims);
      memcpy(dst->data_, data_ + offset, sizeof(Data) * dst->data_size_);
      offset += dst->data_size_;
    }
  }

  void DataStore::map(DataStore* outds, Mapper& m, const View& view)
  {
    check_map_args(outds, view);
    if (data_size_ == 0) {
      return;
    }

    size_t nkeys = 1;
    for (size_t i = 0; i < size_; i++) {
      if (view.dim(i)) {
	nkeys *= value_[i];
      }
    }

    vector< vector<DataPack> > dpgroups(nkeys);

    for (size_t i = 0; i < data_size_; i++) {
      Key tmpkey = index_to_key(i);
      size_t viewed_idx = key_to_viwed_index(tmpkey, view);
      vector<DataPack>& dps = dpgroups.at(viewed_idx);
      dps.push_back(DataPack(tmpkey, &(data_[i])));
    }

    for (size_t i = 0; i < dpgroups.size(); i++) {
      vector<DataPack> &dps = dpgroups.at(i);
      Key viewed_key = key_to_viewed_key(dps.at(0).key(), view);
      m(this, outds, viewed_key, dps);
    }
  }

  void DataStore::load_files(const vector<string>& files, Loader<string>& f)
  {
    load_array(files, f);
  }

  string DataStore::dump(DataPack::Dumper& dumper)
  {
    class WrappedDumper : public Mapper {
    public:
      string result_;
      DataPack::Dumper& dumper_;

      WrappedDumper(DataPack::Dumper& dmpr) : dumper_(dmpr) {}
      int operator()(DataStore *inds, DataStore *outds,
		     Key& key, vector<DataPack>& dps)
      {
	ostringstream os;
	os << "Data Count: " << dps.size() << endl;
	for (vector<DataPack>::iterator itr = dps.begin(); itr != dps.end();
	     itr++) {
	  os << dumper_(*itr);
	}
	result_ = os.str();
	return 0;
      }
    } dmpr(dumper);

    View view(size_);
    for (size_t i = 0; i < size_; i++) {
      view.set_dim(i, false);
    }

    map(NULL, dmpr, view);
    return dmpr.result_;
  }

  void DataStore::set(const size_t *val, Data *dat_ptr)
  {
    if (data_size_ != 0) {
      throw runtime_error("DataStore is already initialized.");
    }

    data_size_ = 1;
    for (size_t i = 0; i < size_; i++) {
      value_[i] = val[i];
      data_size_ *= val[i];
    }
    data_ = dat_ptr;
  }

  size_t DataStore::key_to_index(const Key& key)
  {
    size_t idx = 0;
    for (size_t i = 0; i < size_; i++) {
      size_t offset = 1;
      for (size_t j = i+1; j < size_; j++) {
	offset *= value_[j];
      }
      idx += key.dim(i) * offset;
    }
    return idx;
  }

  Key DataStore::index_to_key(const size_t index)
  {
    Key key(size_);
    size_t _index = index;
    for (size_t i = 0; i < size_; i++) {
      size_t length = 1;
      for (size_t j = i+1; j < size_; j++) {
	length *= value_[j];
      }
      key.set_dim(i, _index / length);
      _index %= length;
    }
    return key;
  }

  size_t DataStore::key_to_viwed_index(const Key& key, const View& view)
  {
    size_t idx = 0;
    for (size_t i = 0; i < size_; i++) {
      if (view.dim(i)) {
	size_t offset = 1;
	for (size_t j = i+1; j < size_; j++) {
	  if (view.dim(j)) {
	    offset *= value_[j];
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
    for (size_t i = 0; i < size_; i++) {
      if (view.dim(i)) {
	viewed_key_size += 1;
      }
    }

    Key viewed_key(viewed_key_size);
    size_t viewed_key_idx = 0;
    for (size_t i = 0; i < size_; i++) {
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
    Key viewed_key(size_);
    for (size_t i = 0; i < size_; i++) {
      if (view.dim(i)) {
	viewed_key.set_dim(i, key.dim(i));
      } else {
	viewed_key.set_dim(i, 0);
      }
    }
    return viewed_key;
  }
#endif

  void DataStore::check_key_range(const Key& key)
  {
    if (size_ != key.size()) {
      throw runtime_error("Dimension size of Key should be same as "
			  "that of DataStore.");
    }
    for (size_t i = 0; i < size_; i++) {
      if (key.dim(i) >= value_[i]) {
	ostringstream os;
	os << "Dimension " << (i+1) << " of Key" << key.to_string()
	   << " is out of range.";
	throw runtime_error(os.str());
      }
    }
  }

  void DataStore::check_map_args(DataStore *outds, const View& view)
  {
    if (this == outds) {
      throw runtime_error("The input and output DataStore should be "
			  "different.");
    }
    if (size_ != view.size()) {
      throw runtime_error("Dimension size of the input DataStore and "
			  "view should be same.");
    }
  }
}
