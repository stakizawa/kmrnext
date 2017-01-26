/// This file is included in data_store.cpp

namespace kmrnext {

  DataPack SimpleFileDataStore::get(const Key& key) {
    load();
    DataPack dp = base::get(key);
    clear_cache();
    return dp;
  }

  vector<DataPack>* SimpleFileDataStore::get(const View& view, const Key& key)
  {
    load();
    vector<DataPack>* dps = base::get(view, key);
    clear_cache();
    return dps;
  }

  DataPack SimpleFileDataStore::remove(const Key& key) {
    load();
    DataPack dp = base::remove(key);
    data_updated_ = true;
    clear_cache();
    return dp;
  }

  string SimpleFileDataStore::filename() {
    ostringstream os;
    os << "./" << this << ".dat";
    return os.str();
  }

}
