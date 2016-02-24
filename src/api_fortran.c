#include "ckmrnext.h"

void *kmrnext_init_ff(void);


void *
kmrnext_init_ff()
{
    printf("this is fortran test.\n");
    return KMRNEXT_init(0, NULL);
}
