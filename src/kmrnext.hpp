#ifndef KMRNEXT_HPP
#define KMRNEXT_HPP
/// \file
/// KMR Next Interface

#include <stdexcept>
#include <sstream>
#include <vector>

#define BACKEND_KMR 1

namespace kmrnext {
  using namespace std;

  /// Maximum dimension size of Key, View and DataStore
  const size_t kMaxDimensionSize = 8;

  class Key;
  class View;
  class Data;
  class DataPack;
  class DataStore;

  ///////////////////////////////////////////////////////////////////////////
  /// A class that stores KMR Next runtime status
  ///////////////////////////////////////////////////////////////////////////
  class KMRNext {
  public:
    /// It initializes the whole system.
    ///
    /// \param[in] argc the number of command line arguments
    /// \param[in] argv the command line arguments
    /// \return         an instance of KMRNext class or its derived class
    static KMRNext* init(int argc, char **argv);

    /// It finalizes the whole system.
    static void finalize();

    /// It creates a DataStore with the specified dimension size.
    ///
    /// \param[in] siz the dimension size of a new DataStore
    /// \return        an instance of DataStore
    DataStore* create_ds(size_t siz);

    KMRNext() {}

    virtual ~KMRNext() {}

  private:
    static KMRNext *kmrnext_;
  };

  ///////////////////////////////////////////////////////////////////////////
  /// Super class of multi-dimensional data class
  ///////////////////////////////////////////////////////////////////////////
  template <typename T>
  class Dimensional {
  protected:
    /// The size of dimension
    size_t size_;
    /// The sizes of each dimension
    T value_[kMaxDimensionSize];

  public:
    /// \param[in] siz the size of dimension
    /// \return        a new instance of Dimensional class
    /// \exception std::runtime_error
    ///                when siz exceeds the maximum dimension size
    explicit Dimensional(size_t siz) : size_(siz) {
      if (size_ > kMaxDimensionSize) {
	ostringstream os;
	os << "Dimension size should be less than " << (kMaxDimensionSize + 1);
	throw runtime_error(os.str());
      }
    }

    virtual ~Dimensional() {};

    virtual void set(const T *val) {
      for (size_t i = 0; i < size_; i++) {
	value_[i] = val[i];
      }
    }

    size_t size() const { return size_; }

    T dim(size_t idx) const {
      if (idx >= size_) {
	ostringstream os;
	os << "Index should be less than dimension size: " << size_;
	throw runtime_error(os.str());
      }
      return value_[idx];
    }

    void set_dim(size_t idx, T val) {
      if (idx >= size_) {
	ostringstream os;
	os << "Index should be less than dimension size: " << size_;
	throw runtime_error(os.str());
      }
      value_[idx] = val;
    }

    string to_string() const {
      ostringstream os;
      os << '<';
      for (size_t i = 0; i < size_; i++) {
	os << value_[i];
	if (i < size_ - 1) {
	  os << ',';
	}
      }
      os << '>';
      return os.str();
    }

    bool operator==(const Dimensional<T>& rhs) const {
      if (size_ != rhs.size_) {
	return false;
      }
      for (size_t i = 0; i < size_; i++) {
	if (value_[i] != rhs.value_[i]) {
	  return false;
	}
      }
      return true;
    }

    bool operator!=(const Dimensional<T>& rhs) const {
      return !(*this == rhs);
    }
  };

  ///////////////////////////////////////////////////////////////////////////
  /// A class that represents key to access data
  ///////////////////////////////////////////////////////////////////////////
  class Key : public Dimensional<size_t> {
  public:
    explicit Key(size_t siz) : Dimensional<size_t>(siz) {}
  };

  ///////////////////////////////////////////////////////////////////////////
  /// A class that defines view of a DataStore
  ///////////////////////////////////////////////////////////////////////////
  class View : public Dimensional<bool> {
  public:
    explicit View(size_t siz) : Dimensional<bool>(siz) {}
  };

  ///////////////////////////////////////////////////////////////////////////
  /// A class that represents a data in a DataStore
  ///////////////////////////////////////////////////////////////////////////
  class Data {
    void *value_;
    size_t value_size_;

  public:
    Data() : value_(NULL), value_size_(0), value_allocated_(false) {}

    Data(void *val, const size_t val_siz)
      : value_(val), value_size_(val_siz), value_allocated_(false) {}

    ~Data();

    void copy_deep(const Data& src);

    /// It returns a pointer to the stored data.
    void *value() { return value_; }

    /// It returns size of the stored data.
    size_t size() { return value_size_; }

  private:
    bool value_allocated_;
  };

  ///////////////////////////////////////////////////////////////////////////
  /// A class that stores a Data ans its Key
  ///////////////////////////////////////////////////////////////////////////
  class DataPack {
  public:
    DataPack(const Key k, Data *d) : key_(k), data_(d) {}

    /// It returns the key.
    Key& key() { return key_; }

    /// It returns the stored data.
    Data *data() { return data_; }

    /////////////////////////////////////////////////////////////////////////
    /// A virtual class used for dumping data in a DataStore
    ///
    /// If you want to dump data in a DataStore, you should implement a class
    ///	that inherites this class and then override operator() method.
    /////////////////////////////////////////////////////////////////////////
    class Dumper {
    public:
      virtual string operator()(DataPack& dp) = 0;
    };

  private:
    Key key_;
    Data *data_;
  };

  ///////////////////////////////////////////////////////////////////////////
  /// A class that stores data
  ///////////////////////////////////////////////////////////////////////////
  class DataStore : public Dimensional<size_t> {
  public:
    explicit DataStore(size_t siz)
      : Dimensional<size_t>(siz), data_(NULL), data_size_(0),
      data_allocated_(false), kmrnext_(NULL) {}

    explicit DataStore(size_t siz, KMRNext *kn)
      : Dimensional<size_t>(siz), data_(NULL), data_size_(0),
      data_allocated_(false), kmrnext_(kn) {}

    virtual ~DataStore();

    /////////////////////////////////////////////////////////////////////////
    /// A virtual class used for applying a function to Data in a DataStore
    ///
    /// If you want to apply a function to each Data in a DataStore (Map),
    /// you should implement a class that inherites this class and then
    /// override operator() method.
    /////////////////////////////////////////////////////////////////////////
    class Mapper {
    public:
      virtual int operator()(DataStore *inds, DataStore *outds, Key& key,
			     vector<DataPack>& dps) = 0;
    };

    /////////////////////////////////////////////////////////////////////////
    /// A virtual class used for loading Data to a DataStore
    ///
    /// If you want to load Data from files or arrays to a DataStore,
    /// you should implement a class that inherites this class and then
    /// override operator() method.
    /////////////////////////////////////////////////////////////////////////
    template <typename Type>
    class Loader {
    public:
      virtual int operator()(DataStore *ds, const Type& param) = 0;
    };

    /// It sets size of each dimension.
    virtual void set(const size_t *val);

    /// It adds a data to this DataStore.
    void add(const Key& key, const Data& data);

    /// It gets a specified data from this DataStore.
    DataPack get(const Key& key);

    /// It gets data whose keys are same when the specified view is applied.
    vector<DataPack>* get(const View& view, const Key& key);

    /// It sets Data from DataStores.
    void set_from(const vector<DataStore*>& dslist);

    /// It splits DataStore to low-dimensional DataStores.
    void split_to(vector<DataStore*>& dslist);

    /// It maps each data.
    void map(DataStore* outds, Mapper& m, const View& view);

    /// It dumps data in the DataStore.
    string dump(DataPack::Dumper& dumper);

    /// It loads files to the DataStore.
    void load_files(const vector<string>& files, Loader<string>& f);

    /// It loads data in an array to the DataStore.
    template <typename Type>
    void load_array(const vector<Type>& array, Loader<Type>& f) {
      if (array.size() == 1) {
	f(this, array[0]);
	return;
      }

      // Check if the size of array is same as the multiple of dimension.
      size_t mi;
      size_t mdim = 1;
      for (mi = 0; mi < size_; mi++) {
	mdim *= value_[mi];
	if (mdim == array.size()) {
	  break;
	}
      }
      if (mdim != array.size()) {
	throw runtime_error("The size of array should be match the "
			    "multiple of dimension sizes of the DataStore.");
      }

      // Calculate size and dimension of Sub DS
      size_t sub_ds_dim = size_ - (mi + 1);
      size_t sub_ds_dims[kMaxDimensionSize];
      size_t sub_ds_siz = 1;
      for (size_t i = 0; i < sub_ds_dim; i++) {
	sub_ds_dims[i] = value_[i + mi + 1];
	sub_ds_siz *= value_[i + mi + 1];
      }

      size_t offset = 0;
      for (size_t i = 0; i < array.size(); i++) {
      	DataStore ds(sub_ds_dim);
      	ds.set(sub_ds_dims, this->data_ + offset);
      	f(&ds, array[i]);
	offset += sub_ds_siz;
      }
    }

  private:
    Data *data_;
    size_t data_size_;
    bool data_allocated_;
    KMRNext *kmrnext_;

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

    // It checks if dimensions of key are inside the range.
    void check_key_range(const Key& key);

    // It checks the arguments of map().
    void check_map_args(DataStore *outds, const View& view);
  };

}

#endif
