#ifndef DATA_STORE_HPP
#define DATA_STORE_HPP
/// \file
/// Basic header for DataStore

#include "kmrnext.hpp"
#include "util.hpp"

#ifdef _OPENMP
#define OMP_FOR_CHUNK_SIZE 4
#endif

namespace kmrnext {
  using namespace std;

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

  /////////////////////////////////////////////////////////////////////////
  // A class that provides utilities for DataStore implementations
  /////////////////////////////////////////////////////////////////////////
  class DataStoreCommonImpl : public DataStore {
    typedef DataStore base;

  public:
    // It returns a Key of the specified indexed Data.
    Key index_to_key(const size_t index);

    // It returns the index of Data calculated from the specified Key.
    size_t key_to_index(const Key& key);

    // It returns the index of Data calculated from the specified Key when
    // the specified view is applied.
    size_t key_to_viewed_index(const Key& key, const View& view);

    // It converts the specified Key by applying the specified View.
    Key key_to_viewed_key(const Key& key, const View& view);

  protected:
    // Index cache
    IndexCache icache_;

    DataStoreCommonImpl(const size_t siz) : base(siz), icache_(IndexCache()) {}

    /// It checks if the given view is correct.
    void check_view(const View& view);

    /// It checks if dimensions of key are inside the range.
    void check_key_range(const Key& key);

    /// It checks the arguments of map().
    void check_map_args(const View& view, DataStore* outds);

    /// It duplicate the DataStore to the specified another one.
    void perform_duplicate(DataStore* ds);

  };

  ///////////////////////////////////////////////////////////////////////////
  /// A template function used for loading array elements into DataStore
  ///////////////////////////////////////////////////////////////////////////
  template <typename T>
  void load_array(const vector<T>& array, DataStore::Loader<T>& loader,
		  KMRNext* next, DataStore* ds,
		  size_t* ds_dims, size_t ds_dims_siz) {
#if VALIDATION
    // Check if the size of array is same as the multiple of dimension.
    {
      size_t prod = 1;
      for (size_t i = 0; i < ds_dims_siz; i++) {
	prod *= ds_dims[i];
	if (prod == array.size()) {
	  break;
	}
      }
      if (array.size() != 1 && array.size() != prod) {
	throw runtime_error("The size of array should be 1 or match the "
			    "product of dimension sizes of the DataStore.");
      }
    }
#endif

    DataStore* ds0 = next->create_ds(1);
    ds0->set_dim(0, array.size());
    Key key(1);
    for (size_t i = 0; i < array.size(); i++) {
      key.set_dim(0, i);
      char* buf;
      size_t buf_siz;
      serialize(array.at(i), &buf, &buf_siz);
      Data dat(buf, buf_siz);
      ds0->add(key, dat);
    }
#ifdef BACKEND_KMR
    // Use Split <T> so that each process loads an array element.
    View split_ds0(1);
    split_ds0.set_dim(0, View::SplitAll);
    ds0->set_split(split_ds0);
#endif

    // Define a mapper for the loader
    class WrappedLoader : public DataStore::Mapper {
    public:
      DataStore::Loader<T>& loader_;

      WrappedLoader(DataStore::Loader<T>& ldr) : loader_(ldr) {}
      int operator()(DataStore* inds, DataStore* outds,
		     Key& k, vector<DataPack>& dps,
		     DataStore::MapEnvironment& env)
      {
	T* val;
	deserialize(static_cast<char*>(dps.at(0).data().value()),
		    dps.at(0).data().size(), &val);
	loader_(outds, *val);
	delete val;
	return 0;
      }
    } wloader(loader);

    View v(1);
    v.set_dim(0, View::SplitAll);
    ds0->map(wloader, v, ds);
    delete ds0;

#ifdef BACKEND_KMR
    View split_ds(ds_dims_siz);
    if (array.size() == 1) {
      // Use Split <All, None, ..> so that data will be distributed to nodes
      // whose count is the top most dimension of the DataStore.
      split_ds.set_dim(0, View::SplitAll);
      for (size_t i = 1; i < ds_dims_siz; i++) {
	split_ds.set_dim(i, View::SplitNone);
      }
    } else {
      // Use Split <All, .., All, None, ..> so that each node stores
      // the contents of each array element.
      long vval = View::SplitAll;
      size_t remain = array.size();
      for (size_t i = 0; i < ds_dims_siz; i++) {
	split_ds.set_dim(i, vval);
	remain /= ds->dim(i);
	if (remain == 1) {
	  vval = View::SplitNone;
	}
      }
    }
    ds->set_split(split_ds);
#endif
  }
}

#endif
