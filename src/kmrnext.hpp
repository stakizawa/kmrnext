#ifndef KMRNEXT_HPP
#define KMRNEXT_HPP
/// \file
/// KMR Next Interface

/// The backend runtime (SERIAL, KMR)
#define BACKEND_KMR 1

#include <stdexcept>
#include <sstream>
#include <vector>
#ifdef BACKEND_KMR
#include <mpi.h>
#include <kmr.h>
#endif

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

    /// It enables profiling option.
    void enable_profile();

    /// It disables profiling option.
    void disable_profile();

    /// It returns true if profiling option is set.
    bool profile() { return profile_; };

#ifdef BACKEND_KMR
    /// It returns MPI processes.
    int nprocs() { return nprocs_; }

    /// It returns rank of this process.
    int rank() { return rank_; }

    /// It returns a KMR object.
    KMR *kmr() { return mr_; };
#endif

  private:
    // True if profiling is on.
    bool profile_;

    static KMRNext *kmrnext_;

    KMRNext();

    ~KMRNext();

#ifdef BACKEND_KMR
    // MPI_Comm where the KMRNext instance belongs
    MPI_Comm world_comm_;
    // Number of processes in the MPI_Comm
    int nprocs_;
    // Rank of this process in the MPI_Comm
    int rank_;
    // KMR instance to be used
    KMR *mr_;
#endif
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
  public:
    Data();

    Data(void *val, const size_t val_siz);

    ~Data();

    void copy_deep(const Data& src);

    void copy_shallow(const Data& src);

    /// It returns a pointer to the stored data.
    void *value() { return value_; }

    /// It returns size of the stored data.
    size_t size() { return value_size_; }

    bool operator==(const Data& rhs) const;

    bool operator!=(const Data& rhs) const {
      return !(*this == rhs);
    }

#ifdef BACKEND_KMR
    /// It returns rank of owner process of this Data.
    int owner() { return owner_; };

    /// It sets owner of this Data.
    void set_owner(int rank) { owner_ = rank; }

    /// It sets that this Data is shared among processes.
    void shared() { shared_ = true; }

    /// It returns true if this Data is shared among processes.
    bool is_shared() { return shared_; }
#endif

  private:
    void *value_;
    size_t value_size_;
    bool value_allocated_;
#ifdef BACKEND_KMR
    int owner_;
    bool shared_;
#endif
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
      data_allocated_(false), parallel_(false),
      kmrnext_(NULL) {}

    explicit DataStore(size_t siz, KMRNext *kn)
      : Dimensional<size_t>(siz), data_(NULL), data_size_(0),
      data_allocated_(false), parallel_(false),
      kmrnext_(kn) {}

    virtual ~DataStore();

    /////////////////////////////////////////////////////////////////////////
    /// A class that stores map function execution environment
    ///
    /// When you call DataStore.Map() to a DataStore, the execution
    /// environment of a Data is set to this class.
    /////////////////////////////////////////////////////////////////////////
    struct MapEnvironment {
      /// Process id.
      ///
      /// In case of Serial backend driver, it is 0.
      /// In case of KMR backend driver, it is a MPI rank in MPI_COMM_WORLD.
      int rank;
#ifdef BACKEND_KMR
      /// MPI_Comm used for processes a Value between processes.
      MPI_Comm mpi_comm;
#endif
    };

    /////////////////////////////////////////////////////////////////////////
    /// A virtual class used for applying a function to Data in a DataStore
    ///
    /// If you want to apply a function to each Data in a DataStore (Map),
    /// you should implement a class that inherites this class and then
    /// override operator() method.
    /////////////////////////////////////////////////////////////////////////
    class Mapper {
    public:
      virtual int operator()(DataStore *inds, DataStore *outds,
			     Key& key, vector<DataPack>& dps,
			     MapEnvironment& env) = 0;
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
    ///
    /// If the data does not exist, it returns NULL data.(TODO comment)
    DataPack get(const Key& key);

    /// It gets data whose keys are same when the specified view is applied.
    vector<DataPack>* get(const View& view, const Key& key);

    /// It removes a specified data from this DataStore.
    ///
    /// If the data does not exist, it returns NULL data.(TODO comment)
    DataPack remove(const Key& key);

    /// It sets Data from DataStores.
    void set_from(const vector<DataStore*>& dslist);

    /// It splits DataStore to low-dimensional DataStores.
    void split_to(vector<DataStore*>& dslist);

    /// It maps each data.
    ///
    /// The mapper object, m, can be run in parallel on data that have
    /// the same key using a given MPI_Comm.
    void map(DataStore* outds, Mapper& m, const View& view);

#ifdef BACKEND_KMR
    /// It maps each data.
    ///
    /// It maps on a single nodes.  Before running the mapper object, m,
    /// data that have the same key are gathered to a node and then the
    /// mapper runs on the nodes as a serial program.
    void map_single(DataStore* outds, Mapper& m, const View& view);

    /// It globally sorts data.
    ///
    /// It changes the arrangement of data elements in the DS among nodes
    /// using the View.  Data elements are distributed among nodes by
    /// dimensions whose values are TURE in the View
    void collate(const View& view);
#endif

    /// It dumps data in the DataStore.
    string dump(DataPack::Dumper& dumper);

    /// It returns stored data count.
    long count();

    /// It loads files to the DataStore.
    ///
    /// \param[in] files the array of files to be loaded
    /// \param[in] f     the function object used to load each files
    /// \exception std::runtime_error
    ///                when there is a mismatch between the number of files
    ///                and dimension sizes of the DataStore.
    void load_files(const vector<string>& files, Loader<string>& f);

    /// It loads data in an array to the DataStore.
    ///
    /// \param[in] array the array of data to be loaded
    /// \param[in] f     the function object used to load each data in
    ///                  the array
    /// \exception std::runtime_error
    ///                when there is a mismatch between the size of array
    ///                and dimension sizes of the DataStore.
    template <typename Type>
    void load_array(const vector<Type>& array, Loader<Type>& f) {
      if (array.size() == 1) {
#ifdef BACKEND_SERIAL
	f(this, array[0]);
#elif defined BACKEND_KMR
	if (kmrnext_->rank() == 0) {
	  f(this, array[0]);
	}
#endif
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

#ifdef BACKEND_SERIAL
      size_t start = 0;
      size_t end = array.size();
      size_t offset = 0;
#elif defined BACKEND_KMR
      // Calculate files assignment to ranks
      size_t quotient = array.size() / kmrnext_->nprocs();
      size_t remain   = array.size() % kmrnext_->nprocs();
      size_t start = kmrnext_->rank() * quotient +
	(((size_t)kmrnext_->rank() < remain)? kmrnext_->rank() : remain);
      size_t end = start + quotient +
	(((size_t)kmrnext_->rank() < remain)? 1 : 0);
      size_t offset = start * sub_ds_siz;
#endif

      for (size_t i = start; i < end; i++) {
      	DataStore ds(sub_ds_dim, this->kmrnext_);
      	ds.set(sub_ds_dims, this->data_ + offset);
	ds.parallel_ = true;
      	f(&ds, array[i]);
	offset += sub_ds_siz;
      }
    }

    /// It returns a Key of the specified indexed Data.
    Key index_to_key(const size_t index);

    /// It converts the specified Key by applying the specified View.
    Key key_to_viewed_key(const Key& key, const View& view);

  private:
    // Pointer to stored Data
    Data *data_;
    // Size of data_
    size_t data_size_;
    // True if the data_ is already allocated
    bool data_allocated_;
    // True if the DataStore should be processed in parallel
    bool parallel_;
    // A KMRNext object that stores execution status
    KMRNext *kmrnext_;

    // It sets size of each dimension.
    // It just sets pointer to data, not performs malloc and memcpy.
    void set(const size_t *val, Data *dat_ptr);

    // It returns the index of Data calculated from the specified Key.
    size_t key_to_index(const Key& key);

    // It returns the index of Data calculated from the specified Key when
    // the specified view is applied.
    size_t key_to_viewed_index(const Key& key, const View& view);

    // It checks if dimensions of key are inside the range.
    void check_key_range(const Key& key);

    // It checks the arguments of map().
    void check_map_args(DataStore *outds, const View& view);

    // It checks the arguments of collate().
    void check_collate_args(const View& view);
  };

}

#endif
