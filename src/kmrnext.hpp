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

    /// It enables profiling option.
    void enable_profile();

    /// It disables profiling option.
    void disable_profile();

    /// It returns true if profiling option is set.
    bool profile() { return profile_; };

    /// It creates a DataStore with the specified dimension size.
    ///
    /// \param[in] siz the dimension size of a new DataStore
    /// \return        an instance of DataStore
    DataStore* create_ds(size_t siz);

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

    /// It deeply copies the specified Data.
    /// If the overwrite option is true, it removes the current data and
    /// overwrites by the specified Data.
    void copy_deep(const Data& src, bool overwrite=false);

    /// It shallowly copies the specified Data.
    void copy_shallow(const Data& src);

    /// It copies the specified buffer to this data.
    void copy_buf(void *val, const size_t val_siz);

    /// It returns a pointer to the stored data.
    void *value() { return value_; }

    /// It returns size of the stored data.
    size_t size() { return value_size_; }

    /// It clears the Data but does not delete it.
    void clear();

    bool operator==(const Data& rhs) const;

    bool operator!=(const Data& rhs) const {
      return !(*this == rhs);
    }

#ifdef BACKEND_KMR
    /// It returns rank of owner process of this Data.
    int owner() { return owner_; }

    /// It sets owner of this Data.
    void set_owner(int rank) { owner_ = rank; }

    /// It sets that this Data is shared among processes.
    void shared() { shared_ = true; }

    /// It sets that this Data is not shared among processes.
    void unshared() { shared_ = false; }

    /// It returns true if this Data is shared among processes.
    bool is_shared() { return shared_; }

    /// It resets the share information.
    void reset_share() { owner_ = -1; shared_ = false; }
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
    DataPack(const Key k, Data *d) : key_(k), data_(d), delete_(false) {}

    ~DataPack() { if (delete_) { delete data_; } }

    /// It returns the key.
    Key& key() { return key_; }

    /// It returns the stored data.
    Data *data() { return data_; }

    /// If it is called, the Data will be deleted when this DataPack
    /// is deleted.
    void set_delete() { delete_ = true; }

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
    bool delete_;
  };

  ///////////////////////////////////////////////////////////////////////////
  /// A class that stores data
  ///////////////////////////////////////////////////////////////////////////
  class DataStore : public Dimensional<size_t> {
  public:
    explicit DataStore(size_t siz);

    explicit DataStore(size_t siz, KMRNext *kn);

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

      /// View used to separate the DataStore
      View view;

#ifdef BACKEND_KMR
      /// Allocation View that define the allocation of data in the DataStore
      View allocation_view;

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
    /// Even if the data does not exist, it returns a DataPack object.
    /// However, the Data of the DataPack is NULL.
    DataPack get(const Key& key);

    /// It gets data whose keys are same when the specified view is applied.
    ///
    /// Even if the data does not exist, it returns a DataPack object.
    /// However, the Data of the DataPack is NULL.
    vector<DataPack>* get(const View& view, const Key& key);

    /// It removes a specified data from this DataStore.
    ///
    /// Even if the data does not exist, it returns a DataPack object.
    /// However, the Data of the DataPack is NULL.
    DataPack remove(const Key& key);

    /// It sets Data from DataStores.
    void set_from(const vector<DataStore*>& dslist);

    /// It splits DataStore to low-dimensional DataStores.
    void split_to(vector<DataStore*>& dslist);

    /// It maps each data.
    ///
    /// The mapper object, m, can be run in parallel on data that have
    /// the same key using a given MPI_Comm.  If the output of the mapper
    /// should be written to another DataStore, specify the last parameter
    /// of the output DataStore, outds.  If the last parameter is omitted,
    /// data elements of this DataStore are updated in-place.
    void map(Mapper& m, const View& view, DataStore* outds=self_);

#ifdef BACKEND_KMR
    /// It sets the allocation view of the DataStore.
    void set_allocation_view(const View& view);

    /// It returns the allocation view of the DataStore.
    View get_allocation_view();

    // It globally sorts data.
    // It changes the arrangement of data elements in the DataStore among
    // nodes using the Allocation View.  It is automatically called in
    // DataStore.map() function, so that explicitly calling this function
    // is not required.
    void collate();
#endif

    /// It dumps data in the DataStore.
    string dump(DataPack::Dumper& dumper);

    /// It returns stored data count.
    long count();

    /// It duplicates the DataStore.
    /// It performs deep copy, that is, all data elements in the DataStore
    /// are deeply copied.  Moreover, it copies all attributes of data
    /// elements.
    DataStore* duplicate();

    /// It loads files to the DataStore.
    /// A file in files is passed to the loader.
    ///
    /// \param[in] files  the array of files to be loaded
    /// \param[in] loader the mapper function object used to load each file
    /// \exception std::runtime_error
    ///                when there is a mismatch between the number of files
    ///                and dimension sizes of the DataStore.
    void load_files(const vector<string>& files, Loader<string>& loader);

    /// It loads integers to the DataStore.
    /// An integer value in ints is passed to the loader.
    ///
    /// \param[in] files  the array of integers to be loaded
    /// \param[in] loader the mapper function object used to load
    ///                   each integer
    /// \exception std::runtime_error
    ///                when there is a mismatch between the number of integers
    ///                and dimension sizes of the DataStore.
    void load_integers(const vector<long>& ints, Loader<long>& loader);

    /// It returns a Key of the specified indexed Data.
    Key index_to_key(const size_t index);

    /// It converts the specified Key by applying the specified View.
    Key key_to_viewed_key(const Key& key, const View& view);

    /// It initializes the static fields.
    /// It should be called before using DataStore class.  If you call
    /// KMRNext::init(), this method is also called.
    static void initialize();

    /// It finalizes the static fields.
    /// It should be called if you no longer use DataStore class.
    /// If you call KMRNext::finalize(), this method is also called.
    static void finalize();

  private:
    // Pointer to stored Data
    Data *data_;
    // Size of data_
    size_t data_size_;
    // True if the data_ is already allocated
    bool data_allocated_;
    // True if the input and output DataStore of map function is same
    bool map_inplace_;
    // True if the DataStore should be processed in parallel
    bool parallel_;
    // A KMRNext object that stores execution status
    KMRNext *kmrnext_;
#ifdef BACKEND_KMR
    /// Allocation View of DataStore, that defines data distribution
    View *allocation_view_;
#endif

    // This is a dummy DataStore that represents this object.
    static DataStore* self_;

    // It sets size of each dimension.
    // It just sets pointer to data, not performs malloc and memcpy.
    void set(const size_t *val, Data *dat_ptr);

    // It returns the index of Data calculated from the specified Key.
    size_t key_to_index(const Key& key);

    // It returns the index of Data calculated from the specified Key when
    // the specified view is applied.
    size_t key_to_viewed_index(const Key& key, const View& view);

    // It checks if the given view is correct.
    void check_view(const View& view);

    // It checks if dimensions of key are inside the range.
    void check_key_range(const Key& key);

    // It checks the arguments of map().
    void check_map_args(const View& view, DataStore* outds);

  };

}

#endif
