/*
  It performs the same operation as test00-basic.cpp.
 */
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include "ckmrnext.h"
#ifdef BACKEND_KMR
#include <mpi.h>
#endif

int rank = 0;

const bool kPrint = true;
const int  kDumpCount = 20;

#define kDimension3 3
const size_t kDim3_0 = 10;
const size_t kDim3_1 = 10;
const size_t kDim3_2 = 10;
#define kDimension2 2
const size_t kDim2_0 = 10;
const size_t kDim2_1 = 10;


void load_data(void *ds);
void print_data_store(void *ds, char *padding, int count);

//////////////////////////////////////////////////////////////////////////////
// Main starts from here.
//////////////////////////////////////////////////////////////////////////////
int
main(int argc, char **argv)
{
    void *next = KMRNEXT_init(argc, argv);
#if BACKEND_KMR
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

    ///////////  Create a DataStore
    void *ds1 = KMRNEXT_create_ds(next, kDimension3);
    size_t sizes3[kDimension3] = {kDim3_0, kDim3_1, kDim3_2};
    KMRNEXT_ds_set_size(ds1, sizes3);
    if (kPrint && rank == 0) {
	printf("0. Create a DataStore\n");
	printf("  DataStore: %s\n", KMRNEXT_ds_string(ds1));
	printf("\n");
    }

    ///////////  Load data contents from a file
    load_data(ds1);
    if (kPrint) {
	if (rank == 0) printf("1. Load data to a DataStore\n");
	print_data_store(ds1, "  ", kDumpCount);
//	print_data_store(ds1, "  ", 0);
	if (rank == 0) printf("\n");
    }

    KMRNEXT_finalize();
    return 0;
}


static int
loader(void *ds, const char *file)
{
    void *key = KMRNEXT_create_key(kDimension3);
    for (size_t i = 0; i < kDim3_0; i++) {
	KMRNEXT_key_set(key, 0, i);
	for (size_t j = 0; j < kDim3_1; j++) {
	    KMRNEXT_key_set(key, 1, j);
	    for (size_t k = 0; k < kDim3_2; k++) {
		KMRNEXT_key_set(key, 2, k);
		long val = (long)(i*j*k);
		void *d = KMRNEXT_create_data(&val, sizeof(long));
		KMRNEXT_ds_add(ds, key, d);
		KMRNEXT_free_data(d);
	    }
	}
    }
    return 0;
}

void
load_data(void *ds)
{
    char *files = {"dummy1"};
    KMRNEXT_ds_load_files(ds, files, 1, loader);
}

static char*
dpdumper(void *dp)
{
    void *key = KMRNEXT_dp_key(dp);
    void *data = KMRNEXT_dp_data(dp);
    const char *key_str = KMRNEXT_key_string(key);
    long *val = (long*)KMRNEXT_data_value(data);
    size_t len = strlen(key_str) + 10;  // "  ", " : ", "000", "\n" and "\0"
    char *buf = (char*)calloc(len, sizeof(char));
    snprintf(buf, len, "  %s : %3ld\n", key_str, *val);
    return buf;  // it may be a memory leak
}

void
print_data_store(void *ds, char *padding, int count)
{
    long ds_count = KMRNEXT_ds_count(ds);
    const char *ds_str = KMRNEXT_ds_dump(ds, dpdumper);
    if (rank != 0) return;

    printf("%sCount of data in the DataStore: %ld\n", padding, ds_count);
    if (count > 0) {
	printf("%sValues (top%d)\n", padding, count);
    } else {
	printf("%sValues (all)\n", padding);
    }

    int cnt = 0;
    char buf[80];
    char *lstart = (char*)ds_str;
    size_t llen = 0;
    for (size_t i = 0; i < strlen(ds_str); i++) {
	if (count > 0 && cnt >= count) {
	    break;
	}
	if (ds_str[i] == '\n') {
	    strncpy(buf, lstart, llen);
	    buf[llen+1] = '\0';
	    printf("%s%s\n", padding, buf);
	    lstart = (char*)&(ds_str[i+1]);
	    llen = 0;
	    cnt += 1;
	} else {
	    llen += 1;
	}
    }
}
