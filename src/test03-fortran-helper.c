#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char* dumper_helper(char *key_str, long val);
void print_data_store_helper(char *dumped_str, long ds_count,
			     int nspace, int print_count);
void print_get_result_helper(char *key_req_str, char *key_ans_str, long val,
			     size_t vsize);
void print_get_view_result_helper1(char *view_str, char*key_str,
				   size_t dp_count, int print_count);
void print_get_view_result_helper2(char *dp_key_str, long dat_val);

#define FORTRAN_PADDING " "


char*
dumper_helper(char *key_str, long val)
{
    // "  ", " : ", "0000", "\n", "\0" and FORTRAN_PADDING
    size_t len = strlen(key_str) + 12;
    char *buf = (char*)calloc(len, sizeof(char));
    snprintf(buf, len, "  %s : %ld\n", key_str, val);
    free(key_str);
    return buf;
}

void
print_data_store_helper(char *dumped_str, long ds_count,
			int nspace, int print_count)
{
    char *padding = (char *)calloc((size_t)nspace + 1, sizeof(char));
    for (int i = 0; i < nspace; i++) {
	padding[i] = ' ';
    }
    padding[nspace] = '\0';

    printf(FORTRAN_PADDING "%sCount of data in the DataStore: %ld\n",
	   padding, ds_count);
    if (print_count > 0) {
	printf(FORTRAN_PADDING "%sValues (top%d)\n", padding, print_count);
    } else {
	printf(FORTRAN_PADDING "%sValues (all)\n", padding);
    }
    int cnt = 0;
    char buf[80];
    char *lstart = (char*)dumped_str;
    size_t llen = 0;
    for (size_t i = 0; i < strlen(dumped_str); i++) {
	if (print_count > 0 && cnt >= print_count) {
	    break;
	}
	if (dumped_str[i] == '\n') {
	    strncpy(buf, lstart, llen);
	    buf[llen] = '\0';
	    printf(FORTRAN_PADDING "%s%s\n", padding, buf);
	    lstart = (char*)&(dumped_str[i+1]);
	    llen = 0;
	    cnt += 1;
	} else {
	    llen += 1;
	}
    }

    free(padding);
    free(dumped_str);
}

void
print_get_result_helper(char *key_req_str, char *key_ans_str, long val,
			size_t vsize)
{
    printf(FORTRAN_PADDING "  Query key: %s    Result: %s: %ld (Size:%ld)\n",
	   key_req_str, key_ans_str, val, vsize);
    free(key_req_str);
    free(key_ans_str);
}

void
print_get_view_result_helper1(char *view_str, char*key_str, size_t dp_count,
			      int print_count)
{
    printf(FORTRAN_PADDING "  Condition\n");
    printf(FORTRAN_PADDING "    view: %s, key: %s\n", view_str, key_str);
    printf(FORTRAN_PADDING "  Result\n");
    printf(FORTRAN_PADDING "    size: %ld\n", dp_count);
    if (print_count > 0) {
	printf(FORTRAN_PADDING "    values (top%d)\n", print_count);
    } else {
	printf(FORTRAN_PADDING "    values (all)\n");
    }
    free(view_str);
    free(key_str);
}

void
print_get_view_result_helper2(char *dp_key_str, long dat_val)
{
    printf(FORTRAN_PADDING "    %s : %ld\n", dp_key_str, dat_val);
    free(dp_key_str);
}
