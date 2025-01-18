#ifndef __UTILS_H__
#define __UTILS_H__

#include "struct_def.h"
#include "lps.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define MARGIN 4

void print_seq(uint64_t id, const char *seq, int seq_len, const char *seq_name, int offset, int rank, FILE *out);

void print_link(uint64_t id1, char sign1, uint64_t id2, char sign2, int overlap, FILE *out);

void print_ref_seq(struct ref_seq *seqs, FILE *out);

int suffix_to_prefix_overlap(const char *str1, const char *str2, int start1, int end1, int start2, int end2);

void free_opt_arg(struct opt_arg *args);

void free_ref_seq(struct ref_seq *seqs);

void find_boundaries(uint64_t start_loc, uint64_t end_loc, struct chr *chrom, uint64_t *latest_end, uint64_t *earliest_start);

void variate(struct chr *chrom, const char *org_seq, const char *alt_token, uint64_t start_loc, int lcp_level, uint64_t *core_id_index, int* failed_var_count, FILE *out, FILE *out_err);

#endif