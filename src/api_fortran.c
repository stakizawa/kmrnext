#include "ckmrnext.h"

typedef struct {
    int rank;
#ifdef BACKEND_KMR
    int mpi_comm;
#endif
    void *p;
} mapenvf;

typedef int (*kmrnext_mapfn_t_ff)(void *ids, void *ods, void *key,
				  datapacks dps, mapenvf env);

void *KMRNEXT_init_ff(void);
void KMRNEXT_ds_map_ff(void *ids, void *ods, void *view,
		       kmrnext_mapfn_t_ff mf, void *p);
#if 0
void print_string(char *str);
void print_strings(char **strs, size_t n);
#endif


void *
KMRNEXT_init_ff()
{
    return KMRNEXT_init(0, NULL);
}

typedef struct {
    kmrnext_mapfn_t_ff m;
    void *p;
} fpackage;

static int
mapfn_wrapper(void *ids, void *ods, void *key, datapacks dps, mapenv env)
{
    fpackage *fp = (fpackage*)env.p;
    mapenvf envf;
    envf.rank = env.rank;
#ifdef BACKEND_KMR
    envf.mpi_comm = MPI_Comm_c2f(env.mpi_comm);
#endif
    envf.p = env.p;
    return (fp->m)(ids, ods, key, dps, envf);
}

void
KMRNEXT_ds_map_ff(void *ids, void *ods, void *view, kmrnext_mapfn_t_ff m,
		  void *p)
{
    fpackage fp;
    fp.m = m;
    fp.p = p;
    KMRNEXT_ds_map(ids, ods, view, mapfn_wrapper, &fp);
}

#if 0
void
print_string(char *str)
{
    printf("Fortran Test: %s\n", str);
}

void
print_strings(char **strs, size_t n)
{
    for (size_t i = 0; i < n; i++) {
	printf("Fortran Test: %s\n", strs[i]);
    }
}
#endif
