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


int summarizer(void *ids, void *ods, void *key, datapacks dps, mapenv env);
void load_data(void *ds);
void print_data_store(void *ds, char *padding, int count);
void print_get_result(void *key, void *dp);
void print_get_view_result(datapacks dps, void *view, void *key, int count);

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
	if (rank == 0) printf("\n");
    }

    ///////////  Setup keys
    void *key1 = KMRNEXT_create_key(kDimension3);
    void *key2 = KMRNEXT_create_key(kDimension3);
    size_t kval1[kDimension3] = {2, 2, 2};
    size_t kval2[kDimension3] = {2, 2, 3};
    KMRNEXT_key_set_size(key1, kval1);
    KMRNEXT_key_set_size(key2, kval2);

    ///////////  Get a data from a DataStore
    void *dp1 = KMRNEXT_ds_get(ds1, key1);
    void *dp2 = KMRNEXT_ds_get(ds1, key2);
    if (kPrint) {
	if (rank == 0) printf("2. Get a data from a DataStore by get()\n");
	print_get_result(key1, dp1);
	print_get_result(key2, dp2);
	if (rank == 0) printf("\n");
    }
    KMRNEXT_free_dp(dp1);
    KMRNEXT_free_dp(dp2);

    ///////////  Setup views
    void *v1 = KMRNEXT_create_view(kDimension3);
    void *v2 = KMRNEXT_create_view(kDimension3);
    void *v3 = KMRNEXT_create_view(kDimension3);
    void *v4 = KMRNEXT_create_view(kDimension3);
    bool flags1[kDimension3] = {true, true, true};
    bool flags2[kDimension3] = {true, false, true};
    bool flags3[kDimension3] = {true, false, false};
    bool flags4[kDimension3] = {false, false, false};
    KMRNEXT_view_set(v1, flags1);
    KMRNEXT_view_set(v2, flags2);
    KMRNEXT_view_set(v3, flags3);
    KMRNEXT_view_set(v4, flags4);

    ///////////  Get a data from a DataStore with a view
    datapacks dps1 = KMRNEXT_ds_get_view(ds1, key1, v1);
    datapacks dps2 = KMRNEXT_ds_get_view(ds1, key2, v1);
    datapacks dps3 = KMRNEXT_ds_get_view(ds1, key1, v2);
    datapacks dps4 = KMRNEXT_ds_get_view(ds1, key2, v2);
    datapacks dps5 = KMRNEXT_ds_get_view(ds1, key1, v3);
    datapacks dps6 = KMRNEXT_ds_get_view(ds1, key2, v3);
    datapacks dps7 = KMRNEXT_ds_get_view(ds1, key1, v4);
    datapacks dps8 = KMRNEXT_ds_get_view(ds1, key2, v4);
    if (kPrint) {
	printf("3. Get data from a DataStore by get(view)\n");
	print_get_view_result(dps1, v1, key1, kDumpCount);
	print_get_view_result(dps2, v1, key2, kDumpCount);
	print_get_view_result(dps3, v2, key1, kDumpCount);
	print_get_view_result(dps4, v2, key2, kDumpCount);
	print_get_view_result(dps5, v3, key1, kDumpCount);
	print_get_view_result(dps6, v3, key2, kDumpCount);
	print_get_view_result(dps7, v4, key1, kDumpCount);
	print_get_view_result(dps8, v4, key2, kDumpCount);
	if (rank == 0) printf("\n");
    }
    KMRNEXT_free_datapacks(dps1);
    KMRNEXT_free_datapacks(dps2);
    KMRNEXT_free_datapacks(dps3);
    KMRNEXT_free_datapacks(dps4);
    KMRNEXT_free_datapacks(dps5);
    KMRNEXT_free_datapacks(dps6);
    KMRNEXT_free_datapacks(dps7);
    KMRNEXT_free_datapacks(dps8);

    ///////////  Apply map functions
    void *ds2 = KMRNEXT_create_ds(next, kDimension2);
    size_t sizes2[kDimension2] = {kDim2_0, kDim2_1};
    KMRNEXT_ds_set_size(ds2, sizes2);
    KMRNEXT_ds_map(ds1, ds2, v2, summarizer);
    if (kPrint) {
	if (rank == 0) printf("4. Apply map to each data in a DataStore\n");
	if (rank == 0) printf("  Output DataStore\n");
	print_data_store(ds2, "    ", kDumpCount);
	if (rank == 0) printf("\n");
    }
    KMRNEXT_free_ds(ds2);

    KMRNEXT_free_view(v1);
    KMRNEXT_free_view(v2);
    KMRNEXT_free_view(v3);
    KMRNEXT_free_view(v4);

    KMRNEXT_free_key(key1);
    KMRNEXT_free_key(key2);

    KMRNEXT_free_ds(ds1);
    KMRNEXT_finalize();
    return 0;
}


int
summarizer(void *ids, void *ods, void *key, datapacks dps, mapenv env)
{
    long sum = 0;
    for (size_t i = 0; i < dps.count; i++) {
	void *dp = dps.data[i];
	long v = *(long*)KMRNEXT_data_value(KMRNEXT_dp_data(dp));
	sum += v;
    }
    void *data = KMRNEXT_create_data(&sum, sizeof(long));
    KMRNEXT_ds_add(ods, key, data);
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
    KMRNEXT_free_key(key);
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
    char *key_str = KMRNEXT_key_string(key);
    long *val = (long*)KMRNEXT_data_value(data);
    size_t len = strlen(key_str) + 11;  // "  ", " : ", "0000", "\n" and "\0"
    char *buf = (char*)calloc(len, sizeof(char));
    snprintf(buf, len, "  %s : %ld\n", key_str, *val);
    free(key_str);
    return buf;  // it may be a memory leak
}

void
print_data_store(void *ds, char *padding, int count)
{
    long ds_count = KMRNEXT_ds_count(ds);
    char *ds_str = KMRNEXT_ds_dump(ds, dpdumper);
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
    free(ds_str);
}

void
print_get_result(void *key, void *dp)
{
    if (rank != 0) return;
    char *key_req_str = KMRNEXT_key_string(key);
    char *key_ans_str = KMRNEXT_key_string(KMRNEXT_dp_key(dp));
    void *data = KMRNEXT_dp_data(dp);
    long val = *(long*)KMRNEXT_data_value(data);
    size_t vsiz = KMRNEXT_data_size(data);
    printf("  Query key: %s    Result: %s: %ld (Size:%ld)\n",
	   key_req_str, key_ans_str, val, vsiz);
    free(key_req_str);
    free(key_ans_str);
}

void
print_get_view_result(datapacks dps, void *view, void *key, int count)
{
    if (rank != 0) return;
    char *view_str = KMRNEXT_view_string(view);
    char *key_str = KMRNEXT_key_string(key);
    printf("  Condition\n");
    printf("    view: %s, key: %s\n", view_str, key_str);
    printf("  Result\n");
    printf("    size: %ld\n", dps.count);
    if (count > 0) {
	printf("    values (top%d)\n", count);
    } else {
	printf("    values (all)\n");
    }
    free(view_str);
    free(key_str);
    int cnt = 0;
    for (size_t i = 0; i < dps.count; i++) {
	if (count > 0 && cnt >= count) {
	    break;
	}
	cnt += 1;
	void *dp = dps.data[i];
	key_str = KMRNEXT_key_string(KMRNEXT_dp_key(dp));
	void *data = KMRNEXT_dp_data(dp);
	long val = *(long*)KMRNEXT_data_value(data);
	printf("    %s : %ld\n", key_str, val);
	free(key_str);
    }
    printf("\n");
}
