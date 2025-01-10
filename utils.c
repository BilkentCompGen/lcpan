#include "utils.h"

void print_seq(uint64_t id, const char *seq, int seq_len, const char *seq_name, int offset, int rank, FILE *out) {
	// printf("S\to%ld\t%s\tSN:Z:%s\tSO:i:%d\tSR:i:%d\n", id, seq, seq_name, offset, rank);
	fprintf(out, "S\to%ld\t", id); 
    fwrite(seq, 1, seq_len, out); 
    fprintf(out, "\tSN:Z:%s\tSO:i:%d\tSR:i:%d\n", seq_name, offset, rank); 
}

void print_link(uint64_t id1, char sign1, uint64_t id2, char sign2, int overlap, FILE *out) {
	fprintf(out, "L\to%ld\t%c\to%ld\t%c\t%dM\n", id1, sign1, id2, sign2, overlap);
}

void print_ref_seq(struct ref_seq *seqs, FILE *out) {

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

	for (int j=0; j<chrom->cores_size; j++) {
		if (cores[j].end < start_loc) {
			*latest_core_index = j;
		}
		if (end_loc < cores[j].start) {
			*first_core_after = j;
			return;
		}
	}
}

void variate(struct chr *chrom, const char *org_seq, const char *alt_token, uint64_t start_loc, int lcp_level, uint64_t *core_id_index, FILE *out) {
	
	// decide on boundaries
	uint64_t end_loc = start_loc + strlen(org_seq);
	uint64_t latest_core_index;
	uint64_t first_core_after;	
	find_boundaries(start_loc, end_loc, chrom, &latest_core_index, &first_core_after);
	latest_core_index = MAX(latest_core_index, 3)-3;
	first_core_after = MIN(first_core_after+3, (uint64_t)chrom->cores_size);

	// get the chromosome
    uint64_t marginal_start = chrom->cores[latest_core_index].start; // start of variation
    uint64_t marginal_end = chrom->cores[first_core_after].end; // end of variation

	// for loop from marginal_start to start_loc the check N
	// for loop from start_ to marginal_end the check N

	// check if there is any N in subseq
	// if there is free(subseg); args->failed++; return;

    uint64_t before_len = start_loc - marginal_start;
    uint64_t alt_len = strlen(alt_token);
    uint64_t after_len = marginal_end - end_loc;
    uint64_t total_len = before_len + alt_len + after_len;

    char *subseq = (char *)malloc((total_len) + 1); // +1 for null terminator
    if (subseq == NULL) {
        fprintf(stderr, "Failed to allocate memory for variated string %ld.\n", total_len);
        return;
    }

    strncpy(subseq, chrom->seq+marginal_start, before_len);

	// printf("start_loc: %ld end_loc: %ld\n", start_loc, end_loc);

	// printf("before: ");
	// fwrite(subseq, 1, before_len, stdout); 
	// printf("\n");

    // Append the alternate sequence
    strcpy(subseq+before_len, alt_token);

	// printf("middle: ");
	// fwrite(subseq+before_len, 1, alt_len, stdout); 
	// printf("\n");

    // Append the 'after variation' sequence
    strncpy(subseq+before_len+alt_len, chrom->seq+end_loc, after_len);

	// printf("after: ");
	// fwrite(subseq+before_len+alt_len, 1, after_len, stdout); 
	// printf("\n");

    // Null-terminate the merged sequence
    subseq[(total_len)] = '\0';

	// printf("final: %s\n", subseq);

	struct lps substr;
	init_lps_offset(&substr, subseq, total_len, marginal_start);
	lps_deepen(&substr, lcp_level);

	// for (int i=0; i<chrom->cores_size; i++) {
	// 	printf("start: %ld, end: %ld, label: %u\n", chrom->cores[i].start, chrom->cores[i].end, chrom->cores[i].label);
	// }
	// printf("\n");

	// for (int i=0; i<substr.size; i++) {
	// 	print_core(&(substr.cores[i]));
	// 	printf(" start: %ld, end: %ld, label: %u", substr.cores[i].start, substr.cores[i].end, substr.cores[i].label);
	// 	printf("\n");
	// }

	const struct simple_core *splitting_point;
	int starting_point = 0; 			// the newly created cores' starting index in lps cores array
	const struct simple_core *merging_core;
	int ending_point = substr.size-1; 	// the newly created cores' ending index in lps cores array

	// refine end
	// we assume here that last core (latest_end and substr.first) it either appear in both
	// or latest_end is not equal to substr.first (SSEQ) and latest_end + 1 is substr.last
	if (chrom->cores[latest_core_index].label != substr.cores[starting_point].label) {
		latest_core_index++;
	}
	if (chrom->cores[latest_core_index].label == substr.cores[starting_point].label) {

		while (chrom->cores[latest_core_index].label == substr.cores[starting_point].label) {
			latest_core_index++;
			starting_point++;
		}

		latest_core_index--;
		splitting_point = &(chrom->cores[latest_core_index]);
		// printf("splitting id o%ld\n", splitting_point->id);
	} else {
		// printf("edge case for start\n");
		latest_core_index--;
		splitting_point = NULL;
	}

	// refine start
	// last logic applies here as well
	if (chrom->cores[first_core_after].label != substr.cores[ending_point].label) {
		first_core_after--;
	}

	if (chrom->cores[first_core_after].label == substr.cores[ending_point].label) {
		
		while (chrom->cores[first_core_after].label == substr.cores[ending_point].label) {
			first_core_after--;
			ending_point--;
		}

		first_core_after++;
		merging_core = &(chrom->cores[first_core_after]);
		// printf("merging id o%ld\n", merging_core->id);
	} else {
		// printf("edge case for end\n");
		first_core_after++;
		merging_core = NULL;
	}

	uint64_t id = *core_id_index;

	if (starting_point <= ending_point) {
		print_seq(id, subseq + substr.cores[starting_point].start - marginal_start, (substr.cores[starting_point].end - substr.cores[starting_point].start), "VAR", substr.cores[starting_point].start, 1, out);
		id++;
	} else {
		// printf("no extra node\n");
	}

	if (splitting_point) {
		print_link(splitting_point->id, '+', id-1, '+', (splitting_point->end - substr.cores[starting_point].start), out);
	}
	
	for (int i=starting_point+1; i<=ending_point; i++) {
		print_seq(id, subseq + substr.cores[i].start - marginal_start, (substr.cores[i].end - substr.cores[i].start), "VAR", substr.cores[i].start, 1, out);
		print_link(id-1, '+', id, '+', (substr.cores[i].end - substr.cores[i].start), out);
		id++;
	}

	if (splitting_point) {
		uint64_t overlap = 0; // todo find overlap at the end
		print_link(id-1, '+', merging_core->id, '+', overlap, out);
	}

	*core_id_index = id;
	
	free_lps(&substr);

	free(subseq);
}