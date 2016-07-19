#include <iostream>
#include <cstdlib>
#include <cstring>
#include "ckmrnext.h"
#include "kmrnext.hpp"

using namespace kmrnext;

char *copy_string_to_cstr(string &str);

void *KMRNEXT_init(int argc, char **argv)
{
  KMRNext *next = KMRNext::init(argc, argv);
  return (void*)next;
}

void KMRNEXT_finalize()
{
  KMRNext::finalize();
}

void KMRNEXT_enable_profile(void *next)
{
  KMRNext *_next = (KMRNext*)next;
  _next->enable_profile();
}

void KMRNEXT_disable_profile(void *next)
{
  KMRNext *_next = (KMRNext*)next;
  _next->disable_profile();
}

bool KMRNEXT_profile(void *next)
{
  KMRNext *_next = (KMRNext*)next;
  return _next->profile();
}

#ifdef BACKEND_KMR
long KMRNEXT_nprocs(void *next)
{
  KMRNext *_next = (KMRNext*)next;
  return (long)_next->nprocs();
}

long KMRNEXT_rank(void *next)
{
  KMRNext *_next = (KMRNext*)next;
  return (long)_next->rank();
}
#endif

void *KMRNEXT_create_ds(void *next, size_t siz)
{
  KMRNext *_next = (KMRNext*)next;
  DataStore *ds = _next->create_ds(siz);
  return (void*)ds;
}

void KMRNEXT_free_ds(void *ds)
{
  DataStore *_ds = (DataStore*)ds;
  delete _ds;
}

void KMRNEXT_ds_set_size(void *ds, size_t *val)
{
  DataStore *_ds = (DataStore*)ds;
  _ds->set(val);
}

void KMRNEXT_ds_load_files(void *ds, char **files, size_t nfiles,
			   kmrnext_loadfn_t l)
{
  vector<string> filevec;
  for (size_t i = 0; i < nfiles; i++) {
    string str(files[i]);
    filevec.push_back(str);
  }
  class WrappedLoader : public DataStore::Loader<string> {
    kmrnext_loadfn_t fn_;
  public:
    WrappedLoader(kmrnext_loadfn_t fn) : fn_(fn) {}
    int operator()(DataStore *_ds, const string& file) {
      const char *file_cstr = file.c_str();
      return fn_((void*)_ds, file_cstr);
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

void *KMRNEXT_ds_get(void *ds, void *key)
{
  DataStore *_ds = (DataStore*)ds;
  Key *_key = (Key*)key;
  DataPack dp = _ds->get(*_key);
  DataPack *_dp = new DataPack(dp.key(), dp.data());
  return (void*)_dp;
}

datapacks KMRNEXT_ds_get_view(void *ds, void *key, void *view)
{
  DataStore *_ds = (DataStore*)ds;
  Key *_key = (Key*)key;
  View *_view = (View*)view;
  vector<DataPack> *dpvec = _ds->get(*_view, *_key);

  datapacks dps;
  dps.count = dpvec->size();
  dps.data = (void**)calloc(dps.count, sizeof(void*));
  size_t idx = 0;
  for (vector<DataPack>::iterator itr = dpvec->begin(); itr != dpvec->end();
       itr++) {
    DataPack *dp = new DataPack((*itr).key(), (*itr).data());
    dps.data[idx] = dp;
    idx += 1;
  }
  delete dpvec;
  return dps;
}

void *KMRNEXT_ds_remove(void *ds, void *key)
{
  DataStore *_ds = (DataStore*)ds;
  Key *_key = (Key*)key;
  DataPack dp = _ds->remove(*_key);
  DataPack *_dp = new DataPack(dp.key(), dp.data());
  return (void*)_dp;
}

void KMRNEXT_ds_map(void *ids, void *ods, void *view, kmrnext_mapfn_t m,
		    void *p)
{
  DataStore *_ids = (DataStore*)ids;
  DataStore *_ods = (DataStore*)ods;
  View *_view = (View*)view;
  class WrappedMapper : public DataStore::Mapper {
    kmrnext_mapfn_t fn_;
    void *p_;
  public:
    WrappedMapper(kmrnext_mapfn_t fn, void *_p) : fn_(fn), p_(_p) {}
    int operator()(DataStore *inds, DataStore *outds, Key& key,
		   vector<DataPack>& dpvec, DataStore::MapEnvironment& env) {
      datapacks dps;
      dps.count = dpvec.size();
      dps.data = (void**)calloc(dps.count, sizeof(void*));
      size_t idx = 0;
      for (vector<DataPack>::iterator itr = dpvec.begin(); itr != dpvec.end();
	   itr++) {
	dps.data[idx] = (void*)&(*itr);
	idx += 1;
      }
      mapenv e;
      e.rank = env.rank;
#ifdef BACKEND_KMR
      e.mpi_comm = env.mpi_comm;
#endif
      e.p = p_;
      int cc = fn_((void*)inds, (void*)outds, (void*)&key, dps, e);
      free(dps.data);
      return cc;
    }
  } mapper(m, p);
  _ids->map(mapper, *_view, _ods);
}

long KMRNEXT_ds_count(void *ds)
{
  DataStore *_ds = (DataStore*)ds;
  return _ds->count();
}

char *KMRNEXT_ds_dump(void *ds, kmrnext_dumpfn_t d)
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
  return copy_string_to_cstr(ds_str);
}

char *KMRNEXT_ds_string(void *ds)
{
  DataStore *_ds = (DataStore*)ds;
  string str = _ds->to_string();
  return copy_string_to_cstr(str);
}

#ifdef BACKEND_KMR
void KMRNEXT_ds_set_split(void *ds, void *split)
{
  DataStore *_ds = (DataStore*)ds;
  View *_split = (View*)split;
  _ds->set_split(*_split);
}

void *KMRNEXT_ds_get_split(void *ds)
{
  DataStore *_ds = (DataStore*)ds;
  View split = _ds->get_split();
  View *_split = new View(split.size());
  for (size_t i = 0; i < split.size(); i++) {
    _split->set_dim(i, split.dim(i));
  }
  return (void*)_split;
}

void KMRNEXT_ds_collate(void *ds)
{
  DataStore *_ds = (DataStore*)ds;
  _ds->collate();
}

bool KMRNEXT_ds_collated(void *ds)
{
  DataStore *_ds = (DataStore*)ds;
  return _ds->collated();
}
#endif

void *KMRNEXT_ds_duplicate(void *ds)
{
  DataStore *_ds = (DataStore*)ds;
  DataStore *_ds2 = _ds->duplicate();
  return (void*)_ds2;
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

void KMRNEXT_key_set_size(void *key, size_t *val)
{
  Key *_key = (Key*)key;
  _key->set(val);
}

void KMRNEXT_key_set(void *key, size_t dim, size_t value)
{
  Key *_key = (Key*)key;
  _key->set_dim(dim, value);
}

char *KMRNEXT_key_string(void *key)
{
  Key *_key = (Key*)key;
  string str = _key->to_string();
  return copy_string_to_cstr(str);
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

size_t KMRNEXT_data_size(void *data)
{
  Data *_data = (Data*)data;
  return _data->size();
}

void *KMRNEXT_create_dp(void *key, void *data)
{
  Key *_key = (Key*)key;
  Data *_data = (Data*)data;
  DataPack *dp = new DataPack(*_key, _data);
  return (void*)dp;
}

void KMRNEXT_free_dp(void *dp)
{
  DataPack *_dp = (DataPack*)dp;
  delete _dp;
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

void *KMRNEXT_create_view(size_t siz)
{
  View *v = new View(siz);
  return (void*)v;
}

void KMRNEXT_free_view(void *view)
{
  View *_view = (View*)view;
  delete _view;
}

void KMRNEXT_view_set(void *view, bool *val)
{
  View *_view = (View*)view;
  _view->set(val);
}

char *KMRNEXT_view_string(void *view)
{
  View *_view = (View*)view;
  string str = _view->to_string();
  return copy_string_to_cstr(str);
}

void KMRNEXT_free_datapacks(datapacks dps)
{
  for (size_t i = 0; i < dps.count; i++) {
    DataPack *dp = (DataPack*)dps.data[i];
    delete dp;
  }
  free(dps.data);
}

char *copy_string_to_cstr(string &str)
{
  size_t len = strlen(str.c_str()) + 1;
  char *buf = new char[len];
  memcpy(buf, str.c_str(), len);
  return buf;
}
