#include <stdio.h>
#include <stdlib.h>
#include <string.h>

char* dumper_helper(char *key_str, long val);
void print_data_store_helper(char *dumped_str, long ds_count,
			     int nspace, int print_count);
void print_get_result_helper(char *key_req_str, char *key_ans_str, long val,
			     size_t vsize);


char*
dumper_helper(char *key_str, long val)
{
    size_t len = strlen(key_str) + 11;  // "  ", " : ", "0000", "\n" and "\0"
    char *buf = (char*)calloc(len, sizeof(char));
    snprintf(buf, len, "  %s : %ld\n", key_str, val);
    return buf; // it may be a memory leak
}

void
print_data_store_helper(char *dumped_str, long ds_count,
			int nspace, int print_count)
{
    char *padding = (char *)calloc((size_t)nspace, sizeof(char));
    for (int i = 0; i < nspace; i++) {
	padding[i] = ' ';
    }

    printf("%sCount of data in the DataStore: %ld\n", padding, ds_count);
    if (print_count > 0) {
	printf("%sValues (top%d)\n", padding, print_count);
    } else {
	printf("%sValues (all)\n", padding);
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
	    printf("%s%s\n", padding, buf);
	    lstart = (char*)&(dumped_str[i+1]);
	    llen = 0;
	    cnt += 1;
	} else {
	    llen += 1;
	}
    }

    free(padding);
}

void
print_get_result_helper(char *key_req_str, char *key_ans_str, long val,
			size_t vsize)
{
    printf("  Query key: %s    Result: %s: %ld (Size:%ld)\n",
	   key_req_str, key_ans_str, val, vsize);
}
