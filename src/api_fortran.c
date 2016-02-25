#include "ckmrnext.h"

void *KMRNEXT_init_ff(void);
/* int KMRNEXT_finalize_ff(void); */
/* int KMRNEXT_free_ds_ff(void *ds); */
/* int KMRNEXT_ds_set_size_ff(void *ds, size_t*val); */

#if 0
void print_string(char *str);
void print_strings(char **strs, size_t n);
#endif


void *
KMRNEXT_init_ff()
{
    printf("C: kmrnext_init_ff.\n"); // TODO
    return KMRNEXT_init(0, NULL);
}

/* int */
/* KMRNEXT_finalize_ff() */
/* { */
/*     printf("C: kmrnext_finalize_ff.\n"); // TODO */
/*     KMRNEXT_finalize(); */
/*     return 0; */
/* } */

/* int */
/* KMRNEXT_free_ds_ff(void *ds) */
/* { */
/*     printf("C: kmrnext_free_ds_ff.\n"); // TODO */
/*     KMRNEXT_free_ds(ds); */
/*     return 0; */
/* } */

/* int */
/* KMRNEXT_ds_set_size_ff(void *ds, size_t *val) */
/* { */
/*     printf("C: kmrnext_ds_set_size_ff.\n"); // TODO */
/*     KMRNEXT_ds_set_size(ds, val); */
/*     return 0; */
/* } */

#if 0
void print_string(char *str)
{
    printf("Fortran Test: %s\n", str);
}

void print_strings(char **strs, size_t n)
{
    for (size_t i = 0; i < n; i++) {
	printf("Fortran Test: %s\n", strs[i]);
    }
}
#endif
