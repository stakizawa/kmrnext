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
  class DataElement;
  class SimpleFileDataStore;
  class SimpleFileDataElement;

  ///////////////////////////////////////////////////////////////////////////
  /// A class that stores KMR Next runtime status
  ///////////////////////////////////////////////////////////////////////////
  class KMRNext {
  public:
    /// IO mode for storing data elements in DataStore
    ///
    /// The default IO mode is Memory.
    enum IOMode {
      Memory = 0,
      File   = 1
    };

    /// It initializes the whole system.
    ///
    /// If KMR backend is used, it also initialize MPI.
    ///
    /// \param[in] argc The number of command line arguments.
    /// \param[in] argv The command line arguments.
    /// \return         An instance of KMRNext class or its derived class.
    static KMRNext* init(int argc, char** argv);

    /// It initializes the whole system.
    ///
    /// \return    An instance of KMRNext class or its derived class.
    static KMRNext* init();

    /// It finalizes the whole system.
    static void finalize();

    /// It aborts the whole system.
    static void abort(int errorcode);

    /// It enables profiling option.
    void enable_profile();

    /// It disables profiling option.
    void disable_profile();

    /// It returns true if profiling option is set.
    bool profile() const { return profile_; };

    /// It sets IO mode of DataSotre.
    ///
    /// \param[in] mode The IO mode.
    void set_io_mode(KMRNext::IOMode mode) { iomode_ = mode; };

    /// It returns the current IO mode.
    KMRNext::IOMode io_mode() const { return iomode_; };

    /// It creates a DataStore with the specified dimension size.
    ///
    /// \param[in] siz The dimension size of a new DataStore
    /// \return        An instance of DataStore
    DataStore* create_ds(size_t siz);

#ifdef BACKEND_KMR
    /// It returns MPI processes.
    int nprocs() const { return nprocs_; }

    /// It returns rank of this process.
    int rank() const { return rank_; }

    /// It returns a KMR object.
    KMR* kmr() const { return mr_; };
#endif

  private:
    // True if profiling is on.
    bool profile_;

    // IO mode of DataStore.
    IOMode iomode_;

    static KMRNext* kmrnext_;

    KMRNext();

    ~KMRNext();

#ifdef BACKEND_KMR
    // MPI_Comm where the KMRNext instance belongs
    MPI_Comm world_comm_;
    // True if MPI is initialized by KMRNext
    bool initiate_mpi;
    // Number of processes in the MPI_Comm
    int nprocs_;
    // Rank of this process in the MPI_Comm
    int rank_;
    // KMR instance to be used
    KMR* mr_;
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
    /// It creates a new dimensional data.
    ///
    /// \param[in] siz The numberof of dimensions
    /// \return        A new instance of Dimensional class
    /// \exception std::runtime_error
    ///                When siz exceeds the maximum dimension size
    explicit Dimensional(const size_t siz) : size_(siz) {
      if (size_ > kMaxDimensionSize) {
	ostringstream os;
	os << "Dimension size should be less than " << (kMaxDimensionSize + 1);
	throw runtime_error(os.str());
      }
    }

    virtual ~Dimensional() {};

    /// It sets value of each dimension.
    ///
    /// \param[in] val An array that stores value of each dimension
    virtual void set(const T* val) {
      for (size_t i = 0; i < size_; i++) {
	value_[i] = val[i];
      }
    }

    /// It returns the number of dimensions.
    size_t size() const { return size_; }

    /// It returns the value of the specified dimension.
    ///
    /// \param[in] idx Index of the target dimension
    /// \return        Value of the dimension
    /// \exception std::runtime_error
    ///                When idx exceeds the maximum dimension size
    T dim(const size_t idx) const {
      if (idx >= size_) {
	ostringstream os;
	os << "Index should be less than dimension size: " << size_;
	throw runtime_error(os.str());
      }
      return value_[idx];
    }

    /// It sets a value to the specified dimension.
    ///
    /// \param[in] idx Index of the target dimension
    /// \param[in] val Value to be set
    /// \exception std::runtime_error
    ///                When idx exceeds the maximum dimension size
    virtual void set_dim(const size_t idx, const T val) {
      if (idx >= size_) {
	ostringstream os;
	os << "Index should be less than dimension size: " << size_;
	throw runtime_error(os.str());
      }
      value_[idx] = val;
    }

    /// It returns a string representation of an instance of this class.
    virtual string to_string() const {
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
    explicit Key(const size_t siz) : Dimensional<size_t>(siz) {}
  };

  ///////////////////////////////////////////////////////////////////////////
  /// A class that defines view of a DataStore
  ///////////////////////////////////////////////////////////////////////////
  class View : public Dimensional<long> {
  public:
    /// A pre-defined constant number that indicates that elements in the
    /// specified dimension are split.
    static const long SplitAll  = -1;

    /// A pre-defined constant number that indicates that elements in the
    /// specified dimension are grouped into one.
    static const long SplitNone = 0;

    explicit View(const size_t siz) : Dimensional<long>(siz) {}

    string to_string() const;
  };

  ///////////////////////////////////////////////////////////////////////////
  /// A class that represents a data in a DataStore
  ///////////////////////////////////////////////////////////////////////////
  class Data {
  public:
    /// It creates a new Data.
    ///
    /// It does not allocate memory space for the value.  If you want to
    /// allocate, call allocate().
    ///
    /// \param[in] val     Value of data to be set.
    /// \param[in] val_siz Size of value.
    Data(void* val, const size_t val_siz)
      : value_(val), value_size_(val_siz), value_allocated_(false) {};

    /// Copy constructor
    ///
    /// The copy constructor do not allocate memory for value.
    Data(const Data& obj);

    ~Data();

    /// It allocates a memory space for the value.
    ///
    /// \exception std::runtime_error
    ///            When an error occurs while allocating.
    void allocate();

    /// It returns a pointer to the value of the Data.
    void* value() const { return value_; }

    /// It returns size of the Data.
    size_t size() const { return value_size_; }

    bool operator==(const Data& rhs) const;

    bool operator!=(const Data& rhs) const {
      return !(*this == rhs);
    }

  private:
    // The actual value of Data.
    void* value_;
    // The size of Data value.
    size_t value_size_;
    // True, if value is allocated.
    bool value_allocated_;
  };

  ///////////////////////////////////////////////////////////////////////////
  /// A class that stores a Data ans its Key
  ///////////////////////////////////////////////////////////////////////////
  class DataPack {
  public:
    /// It creates a new DataPack.
    ///
    /// \param[in] k         Key of the DataPack
    /// \param[in] d         Data of the DataPack
    /// \param[in] allocate  if True, Data is internally allocated
    DataPack(const Key& k, const Data& d, bool allocate=false);

    DataPack(const DataPack &obj);

    /// It returns the key.
    Key& key() { return key_; }

    /// It returns the stored data.
    Data& data() { return data_; }

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
    Data data_;
    bool allocated_;
  };

  ///////////////////////////////////////////////////////////////////////////
  /// A class that stores data
  ///////////////////////////////////////////////////////////////////////////
  class DataStore : public Dimensional<size_t> {
    typedef Dimensional<size_t> base;

  public:
    /// It creates a new DataStore.
    ///
    /// Size of each dimension should be set by DataStore.set() function
    /// before using this DataStore.
    ///
    /// \param[in] size number of dimensions of this DataStore.
    /// \param[in] kn   an instance of KMRNext context
    DataStore(const size_t siz, KMRNext* kn);

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
      /// Split pattern that define the allocation of data in the DataStore
      View split;

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
      virtual int operator()(DataStore* inds, DataStore* outds,
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
      virtual int operator()(DataStore* ds, const Type& param) = 0;
    };

    /// It sets size of each dimension.
    ///
    /// This DataStore is ready to use after calling this function.
    ///
    /// \param[in] val an array that stores value of each dimension
    virtual void set(const size_t* val);

    /// It sets size of a specified dimension of this DataStore.
    ///
    /// \param[in] idx  index of the dimension
    /// \param[in] siz the dimension size
    virtual void set_dim(const size_t idx, const size_t siz);

    /// It sets ZERO data (long interger 0) to all the data elements in
    /// this DataStore.
    void zeroize();

    /// It adds a data to this DataStore.
    ///
    /// \param[in] key  Key in this DataStore where the Data is stored
    /// \param[in] data Data to be added
    /// \exception std::runtime_error when addition failed
    virtual void add(const Key& key, const Data& data);

    /// It gets a specified data from this DataStore.
    ///
    /// \param[in] key Index of Data to be gotten
    /// \return    a DataPack object whose Key is key and value is the gotten
    ///            Data.  As the gotten Data is a copy of Data in the
    ///            DataStore, changing the value of the Data does not take
    ///            any effect.
    ///            Even if the data does not exist, it returns a DataPack
    ///            object.  However, the Data of the DataPack is NULL.
    /// \exception std::runtime_error when failed to get
    virtual DataPack get(const Key& key);

    /// It gets data whose keys are same when the specified view is applied.
    ///
    /// \param[in] view a user specified View of DataStore
    /// \param[in] key  Index of Data to be gotten
    /// \return    a vector of DataPack objects whose Keys are key and values
    ///            are the gotten Data.  As the gotten Data is a copy of Data
    ///            in the DataStore, changing the value of the Data does not
    ///            take any effect.
    ///            Even if the data does not exist, it returns a DataPack
    ///            object.  However, the Data of the DataPack is NULL.
    /// \exception std::runtime_error when failed to get
    virtual vector<DataPack>* get(const View& view, const Key& key);

    /// It removes a specified data from this DataStore.
    ///
    /// \param[in] key Index of Data to be removed
    /// \return    a DataPack object whose Key is key and value is the removed
    ///            Data.  Even if the data does not exist, it returns a
    ///            DataPack object.  However, the Data of the DataPack is NULL.
    /// \exception std::runtime_error when failed to get
    virtual DataPack remove(const Key& key);

    /// It sets Data from DataStores.
    virtual void set_from(const vector<DataStore*>& dslist);

    /// It splits DataStore to low-dimensional DataStores.
    virtual void split_to(vector<DataStore*>& dslist);

    /// It maps each data.
    ///
    /// The mapper object, m, can be run in parallel on data that have
    /// the same key using a given MPI_Comm.  If the output of the mapper
    /// should be written to another DataStore, specify the last parameter
    /// of the output DataStore, outds.  If the last parameter is omitted,
    /// data elements of this DataStore are updated in-place.
    virtual void map(Mapper& m, const View& view, DataStore* outds=self_);

#ifdef BACKEND_KMR
    /// It sets the split pattern of the DataStore.
    void set_split(const View& split);

    /// It returns the split pattern of the DataStore.
    View get_split();

    /// It globally sorts data.
    ///
    /// It changes the arrangement of data elements in the DataStore among
    /// nodes using the Split.  It is automatically called in
    /// DataStore.map() function, so that explicitly calling this function
    /// is not required.
    virtual void collate();

    /// It returns true if the last call of map() or collate() actually
    /// performed collate.  Otherwise it returns false.
    bool collated();
#endif

    /// It dumps data in the DataStore.
    string dump(DataPack::Dumper& dumper);

    /// It returns stored data count.
    long count();

    /// It duplicates the DataStore.
    /// It performs deep copy, that is, all data elements in the DataStore
    /// are deeply copied.  Moreover, it copies all attributes of data
    /// elements.
    virtual DataStore* duplicate();

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

#ifdef BACKEND_KMR
    /// It calles the Loader object on each process in parallel.
    /// It is useful for loading process local data to the DataStore.
    ///
    /// The loader object, loader, called on each process receives
    /// the rank number of the process as the second argument of the
    /// operator() method in long integer.  As any other data is not
    /// passed to the method as arguments, the data to be loaded should
    /// be passed to the loader class as constructor parameters.
    ///
    /// \param[in] loader the loader function object used to load the data
    /// \exception std::runtime_error
    ///                when there is a mismatch between the number of data,
    ///                number of processes and dimension sizes of the
    ///                DataStore.
    void load_parallel(Loader<long>& loader);
#endif

    /// It returns a Key of the specified indexed Data.
    Key index_to_key(const size_t index);

    /// It converts the specified Key by applying the specified View.
    Key key_to_viewed_key(const Key& key, const View& view);

    /// It returns the pointer to the specified indexed DataElement.
    DataElement* data_element_at(const Key& key);

    /// It initializes the static fields.
    ///
    /// It should be called before using DataStore class.  If you call
    /// KMRNext::init(), this method is also called.
    static void initialize(KMRNext* next);

    /// It finalizes the static fields.
    ///
    /// It should be called if you no longer use DataStore class.
    /// If you call KMRNext::finalize(), this method is also called.
    static void finalize();

  protected:
    /////////////////////////////////////////////////////////////////////////
    // A class that caches data for calculating indice and keys.
    /////////////////////////////////////////////////////////////////////////
    class IndexCache {
    public:
      IndexCache();

      // It initializes private members.
      void initialize(const size_t* sizes, const size_t i2k_len,
		      const size_t dim_siz);

      // It returns a key associated with the specified index.
      Key i2k(const size_t index) const;

      // It returns dimension offset of the specified dimension.
      size_t dim_offset(const size_t dim) const;
    private:
      // index to key conversion table
      vector<Key> i2k_table_;
      // dimension offset table
      vector<size_t> doffset_table_;
    };

    // Pointer to stored DataElements
    vector<DataElement*> dlist_;
    // Size of dlist_
    size_t dlist_size_;
    // Index cache
    IndexCache icache_;
    // True if the input and output DataStore of map function is same
    bool map_inplace_;
    // A KMRNext object that stores execution status
    KMRNext* kmrnext_;
#ifdef BACKEND_KMR
    // True if the DataStore should be processed in parallel
    bool parallel_;
    // Split pattern of the DataStore, that defines data distribution
    View* split_;
    // Set to be true if the last call of map() or collate() actually
    // performed collate.
    bool collated_;
#endif

    // This is a dummy DataStore that represents this object.
    static DataStore* self_;

    // It is a real implementation of duplicate().
    void __duplicate(DataStore* ds);

    virtual DataElement* __create_de();

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

#ifdef BACKEND_KMR
    // It returns the index of Data calculated from the specified Key when
    // the specified split is applied.
    size_t key_to_split_index(const Key& key, const View& view);
#endif
  };

  ///////////////////////////////////////////////////////////////////////////
  /// A class that represents an element in a DataStore
  ///////////////////////////////////////////////////////////////////////////
  class DataElement {
  public:
    DataElement();

    virtual ~DataElement();

    /// It sets a Data to the DataElement.
    ///
    /// \param[in] data_value          Pointer to data to be set.
    /// \param[in] data_size           Size of data to be set.
    /// \exception std::runtime_error  When copy failed.
    void set(const void* data_value, const size_t data_size);

    /// It replaces a Data in the DataElement by the specified Data.
    ///
    /// \param[in] data_value          Pointer to data to be set.
    /// \param[in] data_size           Size of data to be set.
    /// \exception std::runtime_error  When copy failed.
    void replace(const void* data_value, const size_t data_size);

    /// It returns data value in the DataElement.
    void* value() { return value_; }

    /// It returns size of data in the DataElement.
    size_t size() { return value_size_; }

    /// It returns true if a Data is set to the DataElement.
    bool is_set() const { return data_set_; }

    /// It clears all attribute of the Data.
    virtual void clear();

#ifdef BACKEND_KMR
    /// It returns rank of owner process of the Data.
    int owner() const { return owner_; }

    /// It sets owner of this Data.
    ///
    /// \param[in] rank The owner rank
    void set_owner(int rank) { owner_ = rank; }

    /// It sets that the Data is shared among processes.
    void shared() { shared_ = true; }

    /// It sets that the Data is not shared among processes.
    void unshared() { shared_ = false; }

    /// It returns true if the Data is shared among processes.
    bool is_shared() const { return shared_; }
#endif

  protected:
    // The actual value of Data.
    char* value_;
    // The size of Data.
    size_t value_size_;
    // True, if Data is set
    bool data_set_;

#ifdef BACKEND_KMR
    // The owner of data in the DataElement.
    int owner_;
    // True if data in the DataElement is shared among processes.
    bool shared_;
#endif

    // It sets a Data to the Dataelement.
    //
    // \param[in] val       Data value to be set
    // \param[in] siz       Size of data to be set
    // \param[in] overwrite If the overwrite option is true, it removes
    //                      the current stored-Data and overwrites it
    //                      by the specified Data.
    // \exception std::runtime_error When copy failed
    virtual void set_data(const void* val, const size_t siz,
			  bool overwrite=false);
  };

  ///////////////////////////////////////////////////////////////////////////
  /// A class that stores data in a file
  ///////////////////////////////////////////////////////////////////////////
  class SimpleFileDataStore : public DataStore {
    typedef DataStore base;

  public:
    /// It creates a new SimpleFileDataStore.
    ///
    /// \param[in] size number of dimensions of this DataStore.
    /// \param[in] kn   an instance of KMRNext context
    SimpleFileDataStore(const size_t siz, KMRNext* kn)
      : base(siz, kn), data_updated_(false), data_cached_(false) {}

    virtual ~SimpleFileDataStore();

    void add(const Key& key, const Data& data);

    DataPack get(const Key& key);

    vector<DataPack>* get(const View& view, const Key& key);

    DataPack remove(const Key& key);

    void set_from(const vector<DataStore*>& dslist);

    void split_to(vector<DataStore*>& dslist);

    void map(Mapper& m, const View& view, DataStore* outds=self_);

    DataStore* duplicate();

#ifdef BACKEND_KMR
    virtual void collate();
#endif

  private:
    // True if the data in DataStore is updated, but not written to a file
    bool data_updated_;
    // True if the Data in a file is cached in memory
    bool data_cached_;

    DataElement* __create_de();

    // It returns name of file that stores data elements of the DataStore.
    string filename();

    // It writes data elements in the DataStore to a file.
    //
    // \return   true if Data are stored in a file
    bool store();

    // It reads data elements in a file and then loads to the DataStore.
    //
    // \return   true if Data are loaded from a file
    bool load();

    // It clears memory cache.
    void clear_cache();
  };

  ///////////////////////////////////////////////////////////////////////////
  /// A class that represents an element in a DataStore that uses files
  ///////////////////////////////////////////////////////////////////////////
  class SimpleFileDataElement : public DataElement {
    typedef DataElement base;

  public:
    SimpleFileDataElement();

    virtual ~SimpleFileDataElement() {}

    /// It clears all attribute of the Data.
    void clear();

    /// It restores the Data in the DataElement from a file buffer.
    ///
    /// \param[in] buf  A file buffer
    void restore(char* buf);

    /// It is called when the Data in the DataElement is written to a file.
    ///
    /// \param[in] start_pos   File start position
    /// \param[in] written_siz Written size in bytes
    void written(size_t start_pos, size_t written_siz);

    /// It clears memory cache of the DataElement.
    void clear_cache();

  private:
    // True, if Data is updated only on memory
    bool data_updated_;
    // The offset of the Data in a file
    size_t data_file_offset_;
    // The size of the Data in a file
    size_t data_file_size_;

    void set_data(const void* val, const size_t siz, bool overwrite=false);
  };

}

#endif
