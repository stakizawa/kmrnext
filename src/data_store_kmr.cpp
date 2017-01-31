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
    size_t vsiz;
    deserialize_collate(kv0.v.p, &vsiz, param->ds, param->rank);
    return MPI_SUCCESS;
  }

  void collate_kmr(KMRNext* kn, DataStore* ds,
		   vector< vector<CollatePack> >& cpgroups, size_t ncpgroups) {
    KMR_KVS* kvs0 = kmr_create_kvs(kn->kmr(),
				   KMR_KV_INTEGER, KMR_KV_OPAQUE);
#if KMR_OMP
#ifdef _OPENMP
#pragma omp parallel for
#endif
#endif
    for (size_t i = 0; i < ncpgroups; i++) {
      vector<CollatePack> &cps = cpgroups.at(i);
      if (cps.size() > 0) {
	for (vector<CollatePack>::iterator itr = cps.begin(); itr != cps.end();
	     itr++) {
	  char* buf = NULL;
	  size_t buf_siz;
	  serialize_collate(&(*itr), &buf, &buf_siz);
	  struct kmr_kv_box kv;
	  kv.klen = static_cast<int>(sizeof(size_t));
	  kv.vlen = static_cast<int>(buf_siz);
	  kv.k.i  = calc_send_target(i, ncpgroups,
				     static_cast<size_t>(kn->nprocs()));
	  kv.v.p  = buf;
	  kmr_add_kv(kvs0, kv);
	  delete[] buf;
	  (*itr).de_->clear();
	}
      }
    }
    kmr_add_kv_done(kvs0);

    KMR_KVS* kvs1 = kmr_create_kvs(kn->kmr(),
				   KMR_KV_INTEGER, KMR_KV_OPAQUE);
    struct kmr_option kmr_shflopt = kmr_noopt;
    kmr_shflopt.key_as_rank = 1;
    kmr_shuffle(kvs0, kvs1, kmr_shflopt);
    //parallel_ = true;
    //map_inplace_ = true;
    param_mapper_collate p0 = { ds, kn->rank() };
    kmr_map(kvs1, NULL, &p0, kmr_noopt, mapper_collate);
    //map_inplace_ = false;
    //parallel_ = false;
  }

  void collate_mpi(KMRNext* kn, DataStore* ds,
		   vector< vector<CollatePack> >& cpgroups, size_t ncpgroups) {
    size_t nprocs = static_cast<size_t>(kn->nprocs());

    // calculate send sizes to each process
    int* sendcounts = new int[nprocs](); // initialized by 0
    vector< vector<CollatePack*> > sendpacks(nprocs);
    // TODO OpenMP parallel
    for (size_t i = 0; i < ncpgroups; i++) {
      vector<CollatePack>& cps = cpgroups.at(i);
      if (cps.size() > 0) {
	size_t dest = calc_send_target(i, ncpgroups, nprocs);
	vector<CollatePack*>& spack = sendpacks.at(dest);
	for (vector<CollatePack>::iterator itr = cps.begin(); itr != cps.end();
	     itr++) {
	  sendcounts[dest] += size_collate(&(*itr));
	  spack.push_back(&(*itr));
	}
      }
    }

    // calculate send displacement of each process
    int* sdispls = new int[nprocs](); // initialized by 0
    size_t total_send_size = sendcounts[0];
    for (size_t i = 1; i < nprocs; i++) {
      sdispls[i] = sdispls[i-1] + sendcounts[i-1];
      total_send_size += sendcounts[i];
    }

    // create send buffer
    char* sendbuf = new char[total_send_size];
    // TODO OpenMP parallel
    for (size_t i = 0; i < nprocs; i++) {
      vector<CollatePack*>& spack = sendpacks.at(i);
      if (spack.size() > 0) {
	size_t offset = 0;
	for (vector<CollatePack*>::iterator itr = spack.begin();
	     itr != spack.end(); itr++) {
	  char* buf = &sendbuf[sdispls[i] + offset];
	  size_t buf_siz;
	  serialize_collate(*itr, &buf, &buf_siz);
	  (*itr)->de_->clear();
	  offset += buf_siz;
	}
      }
    }

    // calculate recv sizes from each process
    int* recvcounts = new int[nprocs](); // initialized by 0
    MPI_Alltoall(sendcounts, 1, MPI_INT, recvcounts, 1, MPI_INT,
		 MPI_COMM_WORLD);

    // calculate recv displacement of each process
    int* rdispls = new int[nprocs](); // initialized by 0
    size_t total_recv_size = recvcounts[0];
    for (size_t i = 1; i < nprocs; i++) {
      rdispls[i] = rdispls[i-1] + recvcounts[i-1];
      total_recv_size += recvcounts[i];
    }

    // create receive buffer
    char* recvbuf = new char[total_recv_size];

    // exchange data
    MPI_Alltoallv(sendbuf, sendcounts, sdispls, MPI_CHAR,
		  recvbuf, recvcounts, rdispls, MPI_CHAR,
		  MPI_COMM_WORLD);

    // load data to DataStore
    // TODO OpenMP parallel
    for (size_t i = 0; i < nprocs; i++) {
      if (recvcounts[i] > 0) {
	for (size_t loaded = 0; loaded < static_cast<size_t>(recvcounts[i]);) {
	  size_t siz_one;
	  deserialize_collate(&recvbuf[rdispls[i] + loaded], &siz_one,
			      ds, kn->rank());
	  loaded += siz_one;
	}
      }
    }

    delete[] sendcounts;
    delete[] sdispls;
    delete[] sendbuf;
    delete[] recvcounts;
    delete[] rdispls;
    delete[] recvbuf;
  }

  int size_collate(const CollatePack* cp) {
    // Memory layout
    // 1. [key]  size    : sizeof(size_t)
    // 2. [key]  content : sizeof(size_t) * count_of_key
    // 3. [data] size    : sizeof(size_t)
    // 4. [data] content : size of 3
    // 5. [data] owner   : sizeof(int)
    //
    // Two parameters of Data, value_allocated_ and shared_ is not serialized.

    size_t size;
    Key k = cp->key_;
    size  = sizeof(size_t);            // 1. [key]  size
    size += sizeof(size_t) * k.size(); // 2. [key]  content
    size += sizeof(size_t);            // 3. [data] size
    size += cp->de_->size();           // 4. [data] content
    size += sizeof(int);               // 5. [data] owner
    return static_cast<int>(size);
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
    if (*buf == NULL) {
      *buf = new char[*buf_siz];
    }
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

  void deserialize_collate(const char* buf, size_t* buf_siz, DataStore* ds,
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
    de.set(p, dat_siz);
    de.set_owner(owner);

    // calculate size
    *buf_siz  = sizeof(size_t);           // 1. [key]  size
    *buf_siz += sizeof(size_t) * key_siz; // 2. [key]  content
    *buf_siz += sizeof(size_t);           // 3. [data] size
    *buf_siz += dat_siz;                  // 4. [data] content
    *buf_siz += sizeof(int);              // 5. [data] owner
  }

}
