#include "utils.h"

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

void print_seq(int tid, uint64_t cid, const char *seq, int seq_len, const char *seq_name, int start, int rank, int is_rgfa, FILE *out) {
	fprintf(out, "S\ts%d.%ld\t", tid, cid);
    fwrite(seq, 1, seq_len, out); 
	if (is_rgfa) {
		fprintf(out, "\tSN:Z:%s\tSO:i:%d\tSR:i:%d\n", seq_name, start, rank); 
	} else {
		fprintf(out, "\n");
	}
}

void print_link(int tid1, uint64_t cid1, char sign1, int tid2, uint64_t cid2, char sign2, uint64_t overlap, int no_overlap, FILE *out) {
    if (!no_overlap) {
        fprintf(out, "L\ts%d.%ld\t%c\ts%d.%ld\t%c\t%ldM\n", tid1, cid1, sign1, tid2, cid2, sign2, overlap);
    } else {
        fprintf(out, "L\ts%d.%ld\t%c\ts%d.%ld\t%c\t%dM\n", tid1, cid1, sign1, tid2, cid2, sign2, 0);
    }
}

void print_ref_seq(const struct ref_seq *seqs, int is_rgfa, int no_overlap, FILE *out) {

	printf("[INFO] Writing reference seq to output file.\n");

	// iterate through each chromosome
	for (int i=0; i<seqs->size; i++) {
		const char *seq_name = seqs->chrs[i].seq_name;
		const char *seq = seqs->chrs[i].seq;
		
		if (seqs->chrs[i].cores_size) {
            const struct simple_core *curr_core = &(seqs->chrs[i].cores[0]);
			uint64_t start = curr_core->start;
			uint64_t end = curr_core->end;
			int seq_len = (int)(end-start);
			print_seq(0, curr_core->id, seq+start, seq_len, seq_name, start, 0, is_rgfa, out);
		}
		
		for (int j=1; j<seqs->chrs[i].cores_size; j++) {
            const struct simple_core *curr_core = &(seqs->chrs[i].cores[j]);
            const struct simple_core *prev_core = &(seqs->chrs[i].cores[j-1]);
			uint64_t curr_start = curr_core->start;
			uint64_t curr_end = curr_core->end;
            uint64_t prev_end = prev_core->end;
			int seq_len = (int)(curr_end-curr_start);
            int overlap = (int)(prev_end-curr_start);
            int offset = 0;
            
            // there might be graps ('N') in genome, hence
            if (prev_end <= curr_start) {
                overlap = 0;
            }

            if (no_overlap) {
                offset = overlap;
            }

            print_seq(0, curr_core->id, seq+curr_start+offset, seq_len-offset, seq_name, curr_start+offset, 0, is_rgfa, out);
            print_link(0, prev_core->id, '+', 0, curr_core->id, '+', overlap, no_overlap, out);
		}
	}

	printf("[INFO] Writing process completed.\n");;
}

uint64_t suffix_to_prefix_overlap(const char *str1, const char *str2, uint64_t start1, uint64_t end1, uint64_t start2, uint64_t end2) {
    uint64_t len_str1 = end1 - start1;
    uint64_t len_str2 = end2 - start2;
    uint64_t max_overlap = 0;

    for (uint64_t i = 1; i <= len_str1 && i <= len_str2; i++) {
        if (strncmp(str1 + end1 - i, str2 + start2, i) == 0) {
            max_overlap = i;
        }
    }
    return max_overlap;
}

void find_boundaries(uint64_t start_loc, uint64_t end_loc, const struct chr *chrom, uint64_t *latest_core_index, uint64_t *first_core_after) {

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

void variate(struct t_arg *t_args, const struct chr *chrom, const char *org_seq, const char *alt_token, uint64_t start_loc) {
	
	// decide on boundaries
	uint64_t end_loc = start_loc+strlen(org_seq);
	uint64_t latest_core_index;
	uint64_t first_core_after;	
	find_boundaries(start_loc, end_loc, chrom, &latest_core_index, &first_core_after);

	if (latest_core_index < MARGIN-1)  {
        t_args->failed_var_count += 1;
        pthread_mutex_lock(t_args->out_err_mutex);
        fprintf(t_args->out_err, "VARIATE-MARGIN-START:\tCHROM: %s,\tPOSITION: %ld,\tORG: %s,\tALT: %s,\tlatest_core_index: %ld\n", chrom->seq_name, start_loc, org_seq, alt_token, latest_core_index);
        fflush(t_args->out_err);
        pthread_mutex_unlock(t_args->out_err_mutex);
		return;
	}
	if (first_core_after+MARGIN-1 >= (uint64_t)chrom->cores_size) {
        t_args->failed_var_count += 1;
        pthread_mutex_lock(t_args->out_err_mutex);
        fprintf(t_args->out_err, "VARIATE-MARGIN-END:\tCHROM: %s,\tPOSITION: %ld,\tORG: %s,\tALT: %s,\tfirst_core_after: %ld,\tcores_size: %d\n", chrom->seq_name, start_loc, org_seq, alt_token, first_core_after, chrom->cores_size);
        fflush(t_args->out_err);
        pthread_mutex_unlock(t_args->out_err_mutex);
		return;
	}

	if (chrom->cores[latest_core_index].end >= start_loc || chrom->cores[latest_core_index+1].end < start_loc) {
        t_args->failed_var_count += 1;
        pthread_mutex_lock(t_args->out_err_mutex);
		fprintf(t_args->out_err, "VARIATE-BOUNDARY:\tIncorrect boundary found.\n");
        fflush(t_args->out_err);
        pthread_mutex_unlock(t_args->out_err_mutex);
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
            t_args->failed_var_count += 1;
            pthread_mutex_lock(t_args->out_err_mutex);
			fprintf(t_args->out_err, "VARIATE-(N at start):\tCHROM: %s,\tPOSITION: %ld,\tORG: %s,\tALT: %s,\tMAR START: %ld,\tMAR END: %ld\n", chrom->seq_name, start_loc, org_seq, alt_token, marginal_start, marginal_end);
            fflush(t_args->out_err);
            pthread_mutex_unlock(t_args->out_err_mutex);
			return;
		}
	}
	// for loop from start_ to marginal_end the check N
	for (uint64_t i = end_loc; i < marginal_end; i++) {
		if (chrom->seq[i] == 'N') {
            t_args->failed_var_count += 1;
            pthread_mutex_lock(t_args->out_err_mutex);
			fprintf(t_args->out_err, "VARIATE-(N at end):\tCHROM: %s,\tPOSITION: %ld,\tORG: %s,\tALT: %s,\tMAR START: %ld,\tMAR END: %ld\n", chrom->seq_name, start_loc, org_seq, alt_token, marginal_start, marginal_end);
			fflush(t_args->out_err);
            pthread_mutex_unlock(t_args->out_err_mutex);
            return;
		}
	}

    uint64_t before_len = start_loc-marginal_start;
    uint64_t alt_len = strlen(alt_token);
    uint64_t after_len = marginal_end-end_loc;
    uint64_t total_len = before_len+alt_len+after_len;

    char *subseq = (char *)malloc(total_len+1); // +1 for null terminator
    if (subseq == NULL) {
        t_args->failed_var_count += 1;
        pthread_mutex_lock(t_args->out_err_mutex);
        fprintf(t_args->out_err, "VARIATE-SUBSEQ:\tCHROM: %s,\tPOSITION: %ld,\tORG: %s,\tALT: %s,\tMAR START: %ld,\tMAR END: %ld\n", chrom->seq_name, start_loc, org_seq, alt_token, marginal_start, marginal_end);
        fflush(t_args->out_err);
        pthread_mutex_unlock(t_args->out_err_mutex);
        return;
    }

    strncpy(subseq, chrom->seq+marginal_start, before_len);
    strcpy(subseq+before_len, alt_token);
    strncpy(subseq+before_len+alt_len, chrom->seq+end_loc, after_len);
    subseq[(total_len)] = '\0';

	struct lps substr;
	init_lps_offset(&substr, subseq, total_len, marginal_start);
	lps_deepen(&substr, t_args->lcp_level);

	// uint64_t debug_latest_core_index = latest_core_index;
	// uint64_t debug_first_core_after = first_core_after;

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

	// if (splitting_core->label != substr.cores[starting_point-1].label ||
	// 	splitting_core->start != substr.cores[starting_point-1].start) {
	// 	fprintf(out_err, "VARIATE-START: Something wrong with starting point\n");
	// 	debug = 1;
	// }
	if (splitting_core->label != substr.cores[starting_point-1].label ||
        splitting_core->start != substr.cores[starting_point-1].start ||
		strncmp(chrom->seq+splitting_core->start, subseq+substr.cores[starting_point-1].start-marginal_start, splitting_core->end-splitting_core->start) != 0) {
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
		debug = 2;
	}

    uint64_t merging_overlap = 0;
	//ending point and merging_core
	uint64_t overlap_start = substr.cores[ending_point].start - marginal_start;
	uint64_t overlap_end = substr.cores[ending_point].end - marginal_start;

	merging_overlap = suffix_to_prefix_overlap(subseq, chrom->seq, overlap_start, overlap_end, merging_core->start, merging_core->end);

    if (merging_overlap == 0) {
		debug = 3;
    }

	if (debug) {
        t_args->failed_var_count += 1;
        pthread_mutex_lock(t_args->out_err_mutex);
        if (debug == 1) {
            // fprintf(t_args->out_err, "VARIATE-START: Chr cores first %d: ", MARGIN);
            // for (uint64_t i=debug_latest_core_index; i<=debug_latest_core_index+3; i++) {
            //     fprintf(t_args->out_err, "core: %u, index: %ld,\t", chrom->cores[i].label, chrom->cores[i].start);
            // }
            // fprintf(t_args->out_err, "Alt Cores first %d: ", MARGIN);
            // for (int i=0; i<4; i++) {
            //     fprintf(t_args->out_err, "core: %u, index: %ld,\t", substr.cores[i].label, substr.cores[i].start);
            // }
            fprintf(t_args->out_err, "VARIATE-START:\tCHROM: %s,\tPOSITION: %ld,\tORG: %s,\tALT: %s,\tMAR START: %ld,\tMAR END: %ld", chrom->seq_name, start_loc, org_seq, alt_token, marginal_start, marginal_end);
        }
        else if (debug == 2) {
            // fprintf(t_args->out_err, "VARIATE-END: Chr cores last %d: ", MARGIN);
            // for (uint64_t i=debug_first_core_after-3; i<=debug_first_core_after; i++) {
            //     fprintf(t_args->out_err, "core: %u, index: %ld,\t", chrom->cores[i].label, chrom->cores[i].start);
            // }
            // fprintf(t_args->out_err, "Alt Cores last %d: ", MARGIN);
            // for (int i=substr.size-4; i<substr.size; i++) {
            //     fprintf(t_args->out_err, "core: %u, index: %ld,\t", substr.cores[i].label, substr.cores[i].start-alt_len+strlen(org_seq));
            // }
            fprintf(t_args->out_err, "VARIATE-END:\tCHROM: %s,\tPOSITION: %ld,\tORG: %s,\tALT: %s,\tMAR START: %ld,\tMAR END: %ld", chrom->seq_name, start_loc, org_seq, alt_token, marginal_start, marginal_end);
        }
        else {
            fprintf(t_args->out_err, "VARIATE-END-OVERLAP:\tCHROM: %s,\tPOSITION: %ld,\tORG: %s,\tALT: %s,\tMAR START: %ld,\tMAR END: %ld", chrom->seq_name, start_loc, org_seq, alt_token, marginal_start, marginal_end);
        }
		fprintf(t_args->out_err, "\n");
		fflush(t_args->out_err);
        pthread_mutex_unlock(t_args->out_err_mutex);
		return;
	}

    t_args->bubble_count += 1;

	// print first splitting node and link from reference graph	
    {
        int splitting_overlap = (splitting_core->end-substr.cores[starting_point].start);
        int seq_start = substr.cores[starting_point].start-marginal_start;
        int seq_len = substr.cores[starting_point].end-substr.cores[starting_point].start;
        int start = substr.cores[starting_point].start;
        int offset = t_args->no_overlap ? splitting_overlap : 0;
        print_seq(t_args->thread_id, t_args->core_id_index, subseq+seq_start+offset, seq_len-offset, "VAR", start+offset, 1, t_args->is_rgfa, t_args->out);
        t_args->core_id_index++;
	    print_link(0, splitting_core->id, '+', t_args->thread_id, t_args->core_id_index-1, '+', splitting_overlap, t_args->no_overlap, t_args->out);
    } 

	for (int i=starting_point+1; i<=ending_point; i++) {
        int overlap = (int)(substr.cores[i-1].end-substr.cores[i].start);
        int seq_start = substr.cores[i].start-marginal_start;
        int seq_len = substr.cores[i].end-substr.cores[i].start;
        int start = substr.cores[i].start;
        int offset = t_args->no_overlap ? overlap : 0;
        print_seq(t_args->thread_id, t_args->core_id_index, subseq+seq_start+offset, seq_len-offset, "VAR", start+offset, 1, t_args->is_rgfa, t_args->out);
        print_link(t_args->thread_id, t_args->core_id_index-1, '+', t_args->thread_id, t_args->core_id_index, '+', overlap, t_args->no_overlap, t_args->out);
		t_args->core_id_index++;
	}

	print_link(t_args->thread_id, t_args->core_id_index-1, '+', 0, merging_core->id, '+', merging_overlap, t_args->no_overlap, t_args->out);
	
	free_lps(&substr);

	free(subseq);

    // if (t_args->core_id_index >> 8 != old_core_id >> 8) {
    //     fflush(t_args->out);
    // }
}