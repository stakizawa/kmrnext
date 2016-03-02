#include "ckmrnext.h"

void *KMRNEXT_init_ff(void);

#if 0
void print_string(char *str);
void print_strings(char **strs, size_t n);
#endif


void *
KMRNEXT_init_ff()
{
    return KMRNEXT_init(0, NULL);
}

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
