#ifndef CKMRNEXT_H
#define CKMRNEXT_H
/* \file KMR Next Interface for C */

/* The backend runtime (SERIAL, KMR) */
#define BACKEND_KMR 1

#include <stdio.h>
#include <stdbool.h>
#ifdef BACKEND_KMR
#include <mpi.h>
#endif

#ifdef __cplusplus
extern "C" {
#endif

const long SPLIT_ALL = -1;
const long SPLIT_NONE = 0;

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
#ifdef BACKEND_KMR
typedef int (*kmrnext_load_parallelfn_t)(void *ds, int rank, void *p);
#endif

void *KMRNEXT_init(int argc, char **argv);
void *KMRNEXT_init0();
void KMRNEXT_finalize();

void KMRNEXT_enable_profile(void *next);
void KMRNEXT_disable_profile(void *next);
bool KMRNEXT_profile(void *next);

#ifdef BACKEND_KMR
long KMRNEXT_nprocs(void *next);
long KMRNEXT_rank(void *next);
#endif

void *KMRNEXT_create_ds(void *next, size_t siz);
void KMRNEXT_free_ds(void *ds);
void KMRNEXT_ds_set_size(void *ds, size_t *val);
void KMRNEXT_ds_zeroize(void *ds);
void KMRNEXT_ds_load_files(void *ds, char **files, size_t nfiles,
			   kmrnext_loadfn_t l);
void KMRNEXT_ds_add(void *ds, void *key, void *data);
void *KMRNEXT_ds_get(void *ds, void *key);
datapacks KMRNEXT_ds_get_view(void *ds, void *key, void *view);
void *KMRNEXT_ds_remove(void *ds, void *key);
void KMRNEXT_ds_map(void *ids, void *ods, void *view, kmrnext_mapfn_t m,
		    void *p);
long KMRNEXT_ds_count(void *ds);
char *KMRNEXT_ds_dump(void *ds, kmrnext_dumpfn_t d);
char *KMRNEXT_ds_string(void *ds);
#ifdef BACKEND_KMR
void KMRNEXT_ds_load_parallel(void *ds, kmrnext_load_parallelfn_t l, void *p);
void KMRNEXT_ds_set_split(void *ds, void *split);
void *KMRNEXT_ds_get_split(void *ds);
void KMRNEXT_ds_collate(void *ds);
bool KMRNEXT_ds_collated(void *ds);
#endif
void *KMRNEXT_ds_duplicate(void *ds);

void *KMRNEXT_create_key(size_t siz);
void KMRNEXT_free_key(void *key);
void KMRNEXT_key_set_size(void *key, size_t *val);
void KMRNEXT_key_set_dim(void *key, size_t dim, size_t value);
size_t KMRNEXT_key_get_dim(void *key, size_t dim);
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
void KMRNEXT_view_set(void *view, long *val);
char *KMRNEXT_view_string(void *view);

void KMRNEXT_free_datapacks(datapacks dps);

#ifdef __cplusplus
}
#endif

#endif
