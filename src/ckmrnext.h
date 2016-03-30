#ifndef CKMRNEXT_H
#define CKMRNEXT_H
/* \file KMR Next Interface for C */

/* The backend runtime (SERIAL, KMR) */
#define BACKEND_SERIAL 1

#include <stdio.h>
#include <stdbool.h>
#ifdef BACKEND_KMR
#include <mpi.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    size_t count;
    void **data;
} datapacks;

typedef struct {
    int rank;
#ifdef BACKEND_KMR
    MPI_Comm mpi_comm;
#endif
    void *p;
} mapenv;

typedef int (*kmrnext_loadfn_t)(void *ds, const char *file);
typedef int (*kmrnext_mapfn_t)(void *ids, void *ods, void *key,
			       datapacks dps, mapenv env);
typedef char* (*kmrnext_dumpfn_t)(void *dp);

void *KMRNEXT_init(int argc, char **argv);
void KMRNEXT_finalize();

void *KMRNEXT_create_ds(void *next, size_t siz);
void KMRNEXT_free_ds(void *ds);
void KMRNEXT_ds_set_size(void *ds, size_t *val);
void KMRNEXT_ds_load_files(void *ds, char **files, size_t nfiles,
			   kmrnext_loadfn_t l);
void KMRNEXT_ds_add(void *ds, void *key, void *data);
void *KMRNEXT_ds_get(void *ds, void *key);
datapacks KMRNEXT_ds_get_view(void *ds, void *key, void *view);
void KMRNEXT_ds_map(void *ids, void *ods, void *view, kmrnext_mapfn_t m,
		    void *p);
long KMRNEXT_ds_count(void *ds);
char *KMRNEXT_ds_dump(void *ds, kmrnext_dumpfn_t d);
char *KMRNEXT_ds_string(void *ds);

void *KMRNEXT_create_key(size_t siz);
void KMRNEXT_free_key(void *key);
void KMRNEXT_key_set_size(void *key, size_t *val);
void KMRNEXT_key_set(void *key, size_t dim, size_t value);
char *KMRNEXT_key_string(void *key);

void *KMRNEXT_create_data(void *val, size_t siz);
void KMRNEXT_free_data(void *data);
void *KMRNEXT_data_value(void *data);
size_t KMRNEXT_data_size(void *data);

void *KMRNEXT_create_dp(void *key, void *data);
void KMRNEXT_free_dp(void *dp);
void *KMRNEXT_dp_key(void *dp);
void *KMRNEXT_dp_data(void *dp);

void *KMRNEXT_create_view(size_t siz);
void KMRNEXT_free_view(void *view);
void KMRNEXT_view_set(void *view, bool *val);
char *KMRNEXT_view_string(void *view);

void KMRNEXT_free_datapacks(datapacks dps);

#ifdef __cplusplus
}
#endif

#endif
