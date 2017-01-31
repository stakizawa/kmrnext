/// This file is inclured in data_store.cpp

namespace kmrnext {

  DataPack SimpleFileDataStore::get(const Key& key) {
#if VALIDATION
    check_key_range(key);
#endif
    load();
    size_t idx = key_to_index(key);
    bool unset = !dlist_[idx].is_set();
    DataPack dp = base::get(key);
    if (unset) {
      // It assume that unset DataElement is set by the above get()
      data_updated_ = true;
      store();
    }
    clear_cache();
    return dp;
  }

  vector<DataPack>* SimpleFileDataStore::get(const View& view, const Key& key)
  {
#if VALIDATION
    check_key_range(key);
#endif
    load();
    vector<DataPack>* dps = base::get(view, key);
    data_updated_ = true; // Though sometimes not updated.
    store();
    clear_cache();
    return dps;
  }

  DataPack SimpleFileDataStore::remove(const Key& key) {
    load();
    DataPack dp0 = base::get(key);
    if (dp0.data().value() != NULL) {
      size_t idx = key_to_index(key);
      dlist_[idx].clear();
      data_updated_ = true;
    }
    clear_cache();
    return dp0;
  }

  void SimpleFileDataStore::collate() {
    // It does not clear cache as callete() is usually called followed by map()
    load();
    base::collate();
    // As data is exchanged between processes, the saved file is stale.
    string fname = filename();
    delete_file(fname);
    store();
  }

  string SimpleFileDataStore::filename() {
    ostringstream os;
    os << "./" << this << "_" << kmrnext_->rank() << ".dat";
    return os.str();
  }

}

namespace kmrtask {

  int mapper_collate(const struct kmr_kv_box kv0, const KMR_KVS* kvi,
		     KMR_KVS* kvo, void* p, const long i) {
    param_mapper_collate* param = static_cast<param_mapper_collate*>(p);
    deserialize_collate(kv0.v.p, kv0.vlen, param->ds, param->rank);
    return MPI_SUCCESS;
  }

  void serialize_collate(const CollatePack* cp, char** buf, size_t* buf_siz) {
    // Memory layout
    // 1. [key]  size    : sizeof(size_t)
    // 2. [key]  content : sizeof(size_t) * count_of_key
    // 3. [data] size    : sizeof(size_t)
    // 4. [data] content : size of 3
    // 5. [data] owner   : sizeof(int)
    //
    // Two parameters of Data, value_allocated_ and shared_ is not serialized.

    // calculate size
    Key k = cp->key_;
    *buf_siz  = sizeof(size_t);            // 1. [key]  size
    *buf_siz += sizeof(size_t) * k.size(); // 2. [key]  content
    *buf_siz += sizeof(size_t);            // 3. [data] size
    *buf_siz += cp->de_->size();           // 4. [data] content
    *buf_siz += sizeof(int);               // 5. [data] owner

    // set data to buf
    *buf = new char[*buf_siz];
    // key
    char* p = *buf;
    *reinterpret_cast<size_t*>(p) = k.size();
    p += sizeof(size_t);
    for (size_t i = 0; i < k.size(); i++) {
      *reinterpret_cast<size_t*>(p) = k.dim(i);
      p += sizeof(size_t);
    }
    // data
    *reinterpret_cast<size_t*>(p) = cp->de_->size();
    p += sizeof(size_t);
    memcpy(p, cp->de_->value(), cp->de_->size());
    p += cp->de_->size();
    *reinterpret_cast<int*>(p) = cp->de_->owner();
  }

  void deserialize_collate(const char* buf, size_t buf_siz, DataStore* ds,
			   int owner) {
    // Memory layout
    // 1. [key]  size    : sizeof(size_t)
    // 2. [key]  content : sizeof(size_t) * count_of_key
    // 3. [data] size    : sizeof(size_t)
    // 4. [data] content : size of 3
    // 5. [data] owner   : sizeof(int)

    // key
    char* p = const_cast<char*>(buf);
    size_t key_siz = *reinterpret_cast<size_t*>(p);
    p += sizeof(size_t);
    size_t* key_ary = new size_t[key_siz];
    for (size_t i = 0; i < key_siz; i++) {
      key_ary[i] = *reinterpret_cast<size_t*>(p);
      p += sizeof(size_t);
    }
    Key k(key_siz);
    k.set(key_ary);
    delete[] key_ary;
    // data
    size_t dat_siz = *reinterpret_cast<size_t*>(p);
    p += sizeof(size_t);
#if 0
    // owner in serialized data is omitted
    p += dat_siz;
    int owner = *reinterpret_cast<int*>(p);
#endif

    DataElement& de = ds->data_element_at(k);
    de.replace(p, dat_siz);
    de.set_owner(owner);
  }

}
