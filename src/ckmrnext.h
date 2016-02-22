#ifndef CKMRNEXT_H
#define CKMRNEXT_H
/* \file KMR Next Interface for C */

/* The backend runtime (SERIAL, KMR) */
#define BACKEND_SERIAL 1

#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef int (*kmrnext_loadfn_t)(void *ds, const char *file);
typedef char* (*kmrnext_dumpfn_t)(void *dp);

void *KMRNEXT_init(int argc, char **argv);
void KMRNEXT_finalize();

void *KMRNEXT_create_ds(void *next, size_t siz);
void KMRNEXT_ds_set_size(void *ds, size_t *val);
void KMRNEXT_ds_load_files(void *ds, char *files, size_t nfiles,
			   kmrnext_loadfn_t l);
void KMRNEXT_ds_add(void *ds, void *key, void *data);
long KMRNEXT_ds_count(void *ds);
const char *KMRNEXT_ds_dump(void *ds, kmrnext_dumpfn_t d);
const char *KMRNEXT_ds_string(void *ds);

void *KMRNEXT_create_key(size_t siz);
void KMRNEXT_free_key(void *key);
void KMRNEXT_key_set(void *key, size_t dim, size_t value);
const char *KMRNEXT_key_string(void *key);

void *KMRNEXT_create_data(void *val, size_t siz);
void KMRNEXT_free_data(void *data);
void *KMRNEXT_data_value(void *data);

void *KMRNEXT_dp_key(void *dp);
void *KMRNEXT_dp_data(void *dp);

#ifdef __cplusplus
}
#endif

#endif
