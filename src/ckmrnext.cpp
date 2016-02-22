#include "ckmrnext.h"
#include "kmrnext.hpp"
#include <iostream>

using namespace kmrnext;

void *KMRNEXT_init(int argc, char **argv)
{
  KMRNext *next = KMRNext::init(argc, argv);
  return (void*)next;
}

void KMRNEXT_finalize()
{
  KMRNext::finalize();
}

void *KMRNEXT_create_ds(void *next, size_t siz)
{
  KMRNext *_next = (KMRNext*)next;
  DataStore *ds = _next->create_ds(siz);
  return (void*)ds;
}

void KMRNEXT_ds_set_size(void *ds, size_t *val)
{
  DataStore *_ds = (DataStore*)ds;
  _ds->set(val);
}

void KMRNEXT_ds_load_files(void *ds, char *files, size_t nfiles,
			   kmrnext_loadfn_t l)
{
  vector<string> filevec;
  for (size_t i = 0; i < nfiles; i++) {
    string str(&files[i]);
    filevec.push_back(str);
  }
  class WrappedLoader : public DataStore::Loader<string> {
    kmrnext_loadfn_t fn_;
  public:
    WrappedLoader(kmrnext_loadfn_t fn) : fn_(fn) {}
    int operator()(DataStore *ds, const string& file) {
      const char *file_cstr = file.c_str();
      return fn_((void*)ds, file_cstr);
    }
  } loader(l);
  DataStore *_ds = (DataStore*)ds;
  _ds->load_files(filevec, loader);
}

void KMRNEXT_ds_add(void *ds, void *key, void *data)
{
  DataStore *_ds = (DataStore*)ds;
  Key *_key = (Key*)key;
  Data *_data = (Data*)data;
  _ds->add(*_key, *_data);
}

long KMRNEXT_ds_count(void *ds)
{
  DataStore *_ds = (DataStore*)ds;
  return _ds->count();
}

const char *KMRNEXT_ds_dump(void *ds, kmrnext_dumpfn_t d)
{
  class WrappedDumper : public DataPack::Dumper {
    kmrnext_dumpfn_t fn_;
  public:
    WrappedDumper(kmrnext_dumpfn_t fn) : fn_(fn) {}
    string operator()(DataPack& dp) {
      char * str_c = fn_((void*)&dp);
      return string(str_c);
    }
  } dumper(d);
  DataStore *_ds = (DataStore*)ds;
  string ds_str = _ds->dump(dumper);
  return ds_str.c_str();
}

const char *KMRNEXT_ds_string(void *ds)
{
  DataStore *_ds = (DataStore*)ds;
  string str = _ds->to_string();
  return str.c_str();
}

void *KMRNEXT_create_key(size_t siz)
{
  Key *key = new Key(siz);
  return (void*)key;
}

void KMRNEXT_free_key(void *key)
{
  Key *_key = (Key*)key;
  delete _key;
}

void KMRNEXT_key_set(void *key, size_t dim, size_t value)
{
  Key *_key = (Key*)key;
  _key->set_dim(dim, value);
}

const char *KMRNEXT_key_string(void *key)
{
  Key *_key = (Key*)key;
  string str = _key->to_string();
  return str.c_str();
}

void *KMRNEXT_create_data(void *val, size_t siz)
{
  Data *d = new Data(val, siz);
  return (void*)d;
}

void KMRNEXT_free_data(void *data)
{
  Data *_data = (Data*)data;
  delete _data;
}

void *KMRNEXT_data_value(void *data)
{
  Data *_data = (Data*)data;
  return _data->value();
}

void *KMRNEXT_dp_key(void *dp)
{
  DataPack *_dp = (DataPack*)dp;
  return (void*)(&_dp->key());
}

void *KMRNEXT_dp_data(void *dp)
{
  DataPack *_dp = (DataPack*)dp;
  Data *data = _dp->data();
  return (void*)data;
}
