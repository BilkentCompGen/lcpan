#include "utils.h"

void print_seq(uint64_t id, const char *seq, int seq_len, const char *seq_name, int offset, int rank, FILE *out) {
	fprintf(out, "S\ts%ld\t", id); 
    fwrite(seq, 1, seq_len, out); 
    fprintf(out, "\tSN:Z:%s\tSO:i:%d\tSR:i:%d\n", seq_name, offset, rank); 
}

void print_link(uint64_t id1, char sign1, uint64_t id2, char sign2, int overlap, FILE *out) {
	fprintf(out, "L\ts%ld\t%c\ts%ld\t%c\t%dM\n", id1, sign1, id2, sign2, overlap);
}

void print_ref_seq(struct ref_seq *seqs, FILE *out) {

	printf("[INFO] Writing reference seq to output file.\n");

	// iterate through each chromosome
	for (int i=0; i<seqs->size; i++) {
		const char *seq_name = seqs->chrs[i].seq_name;
		const char *seq = seqs->chrs[i].seq;
		
		if (seqs->chrs[i].cores_size) {
			uint64_t start = seqs->chrs[i].cores[0].start;
			uint64_t end = seqs->chrs[i].cores[0].end;
			int seq_len = (int)(end-start);
			print_seq(seqs->chrs[i].cores[0].id, seq+start, seq_len, seq_name, start, 0, out);
		}
		
		for (int j=1; j<seqs->chrs[i].cores_size; j++) {
			uint64_t start = seqs->chrs[i].cores[j].start;
			uint64_t end = seqs->chrs[i].cores[j].end;
			int seq_len = (int)(end-start);
			
			print_seq(seqs->chrs[i].cores[j].id, seq+start, seq_len, seq_name, start, 0, out);
			print_link(seqs->chrs[i].cores[j-1].id, '+', seqs->chrs[i].cores[j].id, '+', (int)(seqs->chrs[i].cores[j-1].end-seqs->chrs[i].cores[j].start), out);
		}
	}

	printf("[INFO] Writing process completed.\n");
}

int suffix_to_prefix_overlap(const char *str1, const char *str2, int start1, int end1, int start2, int end2) {
    int len_str1 = end1 - start1;
    int len_str2 = end2 - start2;
    int max_overlap = 0;

    for (int i = 1; i <= len_str1 && i <= len_str2; i++) {
        if (strncmp(str1 + end1 - i, str2 + start2, i) == 0) {
            max_overlap = i;
        }
    }
    return max_overlap;
}

void free_opt_arg(struct opt_arg *args) {
    free(args->fasta_fai_path);
	// the rest of the args char * will be freed by getops. hence, no need to free them
}

void free_ref_seq(struct ref_seq *seqs) {
	if (seqs->size) {
		for (int i=0; i<seqs->size; i++) {
			free(seqs->chrs[i].seq_name);
			free(seqs->chrs[i].seq);
			free(seqs->chrs[i].cores);
		}
		free(seqs->chrs);
		seqs->size = 0;
	}
}

void find_boundaries(uint64_t start_loc, uint64_t end_loc, struct chr *chrom, uint64_t *latest_core_index, uint64_t *first_core_after) {

	*latest_core_index = 0;
	*first_core_after = UINT64_MAX;
	const struct simple_core *cores = chrom->cores;	
	int left, mid, right;

	left = 0;
    right = chrom->cores_size;
    while (left < right) {
        mid = left + (right - left) / 2;
        if (cores[mid].end < start_loc) {
            *latest_core_index = mid;
            left = mid + 1;
        } else {
            right = mid;
        }
    }

	*first_core_after = *latest_core_index;
	for (int j=*latest_core_index; j<chrom->cores_size; j++) {
		if (end_loc < cores[j].start) {
			*first_core_after = j;
			return;
		}
	}
}

void variate(struct chr *chrom, const char *org_seq, const char *alt_token, uint64_t start_loc, int lcp_level, uint64_t *core_id_index, int* failed_var_count, FILE *out, FILE *out_err) {
	
	// decide on boundaries
	uint64_t end_loc = start_loc+strlen(org_seq);
	uint64_t latest_core_index;
	uint64_t first_core_after;	
	find_boundaries(start_loc, end_loc, chrom, &latest_core_index, &first_core_after);

	if (latest_core_index < MARGIN-1)  {
		*failed_var_count = *failed_var_count + 1;
		fprintf(out_err, "VARIATE-MARGIN: %s\t%ld\t%s\t%s\n", chrom->seq_name, start_loc, org_seq, alt_token);
		return;
	}
	if (first_core_after+MARGIN-1 >= (uint64_t)chrom->cores_size) {
		*failed_var_count = *failed_var_count + 1;
		fprintf(out_err, "VARIATE-MARGIN: %s\t%ld\t%s\t%s\n", chrom->seq_name, start_loc, org_seq, alt_token);
		return;
	}

	if (chrom->cores[latest_core_index].end >= start_loc || chrom->cores[latest_core_index+1].end < start_loc) {
		*failed_var_count = *failed_var_count + 1;
		fprintf(out_err, "VARIATE-BOUNDARY: Incorrect boundary found.\n");
		return;
	}

	latest_core_index -= MARGIN-1;
	first_core_after += MARGIN-1;

	// get the chromosome
    uint64_t marginal_start = chrom->cores[latest_core_index].start; // start of variation
    uint64_t marginal_end = chrom->cores[first_core_after].end; // end of variation

	// for loop from marginal_start to start_loc the check N
	for (uint64_t i = marginal_start; i < start_loc; i++) {
		if (chrom->seq[i] == 'N') {
			*failed_var_count = *failed_var_count + 1;
			fprintf(out_err, "%s\t%ld\t%s\t%s\n", chrom->seq_name, start_loc, org_seq, alt_token);
			return;
		}
	}
	// for loop from start_ to marginal_end the check N
	for (uint64_t i = end_loc; i < marginal_end; i++) {
		if (chrom->seq[i] == 'N') {
			*failed_var_count = *failed_var_count + 1;
			fprintf(out_err, "%s\t%ld\t%s\t%s\n", chrom->seq_name, start_loc, org_seq, alt_token);
			return;
		}
	}

    uint64_t before_len = start_loc-marginal_start;
    uint64_t alt_len = strlen(alt_token);
    uint64_t after_len = marginal_end-end_loc;
    uint64_t total_len = before_len+alt_len+after_len;

    char *subseq = (char *)malloc(total_len+1); // +1 for null terminator
    if (subseq == NULL) {
		*failed_var_count = *failed_var_count + 1;
        fprintf(out_err, "VARIATE-SUBSEQ: Failed to allocate memory for variated string %ld.\n", total_len);
        return;
    }

    strncpy(subseq, chrom->seq+marginal_start, before_len);
    strcpy(subseq+before_len, alt_token);
    strncpy(subseq+before_len+alt_len, chrom->seq+end_loc, after_len);
    subseq[(total_len)] = '\0';

	struct lps substr;
	init_lps_offset(&substr, subseq, total_len, marginal_start);
	lps_deepen(&substr, lcp_level);

	uint64_t debug_latest_core_index = latest_core_index;
	uint64_t debug_first_core_after = first_core_after;

	int debug = 0;

	// refine start
	// we are using the locations (start) of the cores for the refinements,
	// and matching the ones that are identified in both alt and org sequences.
	int starting_point = 0; 			// the newly created cores' starting index in lps cores array
	while (chrom->cores[latest_core_index].end < start_loc) {
		if (chrom->cores[latest_core_index].start == substr.cores[starting_point].start &&
			chrom->cores[latest_core_index].label == substr.cores[starting_point].label) {
			break; // found first match
		} else if (chrom->cores[latest_core_index].start < substr.cores[starting_point].start) {
			latest_core_index++;
		} else {
			starting_point++;
		}
	}
	// go until there is a mismatch
	while ( chrom->cores[latest_core_index].end<start_loc && 
			chrom->cores[latest_core_index].label == substr.cores[starting_point].label &&
			chrom->cores[latest_core_index].start == substr.cores[starting_point].start) {
		latest_core_index++;
		starting_point++;
	}
	latest_core_index--; // as there should be at least one match
	const struct simple_core *splitting_core = &(chrom->cores[latest_core_index]);

	if (splitting_core->label != substr.cores[starting_point-1].label ||
		splitting_core->start != substr.cores[starting_point-1].start) {
		fprintf(out_err, "VARIATE-START: Something wrong with starting point\n");
		debug = 1;
	}

	// refine end
	int ending_point = substr.size-1; 	// the newly created cores' ending index in lps cores array
	while (end_loc <= chrom->cores[first_core_after].start) {
		if (chrom->cores[first_core_after].label == substr.cores[ending_point].label &&
			chrom->cores[first_core_after].start == substr.cores[ending_point].start-alt_len+strlen(org_seq)) {
			break; // found first match
		} else if (chrom->cores[first_core_after].start < substr.cores[ending_point].start-alt_len+strlen(org_seq)) {
			ending_point--;
		} else {
			first_core_after--;
		}
	}
	// go until there is a mismatch
	while ( end_loc <= chrom->cores[first_core_after].start && 
			chrom->cores[first_core_after].label == substr.cores[ending_point].label && 
			chrom->cores[first_core_after].start == substr.cores[ending_point].start-alt_len+strlen(org_seq)) {
		first_core_after--;
		ending_point--;
	}
	first_core_after++;
	const struct simple_core *merging_core = &(chrom->cores[first_core_after]);

	if (merging_core->label != substr.cores[ending_point+1].label || 
		merging_core->start != substr.cores[ending_point+1].start-alt_len+strlen(org_seq)) {
		fprintf(out_err, "VARIATE-END: Something wrong with ending point\n");
		debug = 1;
	}

	if (debug) {
		fprintf(out_err, "Chr cores first %d: \n", MARGIN);
		for (uint64_t i=debug_latest_core_index; i<=debug_latest_core_index+3; i++) {
			fprintf(out_err, "core: %u, index: %ld\n", chrom->cores[i].label, chrom->cores[i].start);
		}
		fprintf(out_err, "Alt Cores first %d: \n", MARGIN);
		for (int i=0; i<4; i++) {
			fprintf(out_err, "core: %u, index: %ld\n", substr.cores[i].label, substr.cores[i].start);
		}

		fprintf(out_err, "Chr cores last %d: \n", MARGIN);
		for (uint64_t i=debug_first_core_after-3; i<=debug_first_core_after; i++) {
			fprintf(out_err, "core: %u, index: %ld\n", chrom->cores[i].label, chrom->cores[i].start);
		}
		fprintf(out_err, "Alt Cores last %d: \n", MARGIN);
		for (int i=substr.size-4; i<substr.size; i++) {
			fprintf(out_err, "core: %u, index: %ld\n", substr.cores[i].label, substr.cores[i].start-alt_len+strlen(org_seq));
		}

		fprintf(out_err, "\n");
		
		*failed_var_count = *failed_var_count + 1;
		return;
	}

	uint64_t id = *core_id_index;

	// print first splitting node and link from reference graph
	print_seq(id, subseq+substr.cores[starting_point].start-marginal_start, substr.cores[starting_point].end-substr.cores[starting_point].start, "VAR", substr.cores[starting_point].start, 1, out);
	id++;
	print_link(splitting_core->id, '+', id-1, '+', (splitting_core->end-substr.cores[starting_point].start), out);
	
	for (int i=starting_point+1; i<=ending_point; i++) {
		print_seq(id, subseq + substr.cores[i].start - marginal_start, substr.cores[i].end - substr.cores[i].start, "VAR", substr.cores[i].start, 1, out);
		print_link(id-1, '+', id, '+', (substr.cores[i-1].end-substr.cores[i].start), out);
		id++;
	}

	uint64_t overlap = 0;
	//ending point and merging_core
	int overlap_start = substr.cores[ending_point].start - marginal_start;
	int overlap_end = substr.cores[ending_point].end - marginal_start;

	overlap = suffix_to_prefix_overlap(subseq, chrom->seq, overlap_start, overlap_end, merging_core->start, merging_core->end);

	print_link(id-1, '+', merging_core->id, '+', overlap, out);

	*core_id_index = id;
	
	free_lps(&substr);

	free(subseq);
}