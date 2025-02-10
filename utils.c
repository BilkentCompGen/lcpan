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
    if (no_overlap) {
        fprintf(out, "L\ts%d.%ld\t%c\ts%d.%ld\t%c\t%dM\n", tid1, cid1, sign1, tid2, cid2, sign2, 0);
    } else {
        fprintf(out, "L\ts%d.%ld\t%c\ts%d.%ld\t%c\t%ldM\n", tid1, cid1, sign1, tid2, cid2, sign2, overlap);
    }
}

void print_ref_seq(const struct ref_seq *seqs, int is_rgfa, int no_overlap, FILE *out) {

	printf("[INFO] Writing reference seq to output file.\n");

    fprintf(out, "H\tVN:Z:1.1\n");
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

	printf("[INFO] Writing process completed.\n");
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
    *first_core_after = chrom->cores_size;
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

void variate_snp(struct t_arg *t_args, const struct chr *chrom, const char *alt_token, uint64_t start_loc, uint64_t end_loc, uint64_t latest_core_index, uint64_t first_core_after, uint64_t marginal_start, uint64_t marginal_end, uint64_t end_overlap) {

	// print new node connecting directly from latest core before and first core after the altenrating token
    // print new node
    fprintf(t_args->out, "S\ts%d.%ld\t", t_args->thread_id, t_args->core_id_index);
    fwrite(chrom->seq+marginal_start, 1, start_loc-marginal_start, t_args->out);
    fprintf(t_args->out, "%s", alt_token);
    fwrite(chrom->seq+end_loc, 1, marginal_end-end_loc, t_args->out); 
	if (t_args->is_rgfa) {
		fprintf(t_args->out, "\tSN:Z:%d.%ld\tSO:i:%ld\tSR:i:1\n", t_args->thread_id, t_args->core_id_index, marginal_start); 
	} else {
		fprintf(t_args->out, "\n");
	}
    // print splitting link
    print_link(0, chrom->cores[latest_core_index].id, '+', t_args->thread_id, t_args->core_id_index, '+', 0, t_args->no_overlap, t_args->out);
	// print merging link
    print_link(t_args->thread_id, t_args->core_id_index, '+', 0, chrom->cores[first_core_after].id, '+', end_overlap, t_args->no_overlap, t_args->out);
    t_args->core_id_index++;
}

void variate_sv(struct t_arg *t_args, const struct chr *chrom, const char *alt_token, uint64_t start_loc, uint64_t end_loc, uint64_t latest_core_index, uint64_t first_core_after, uint64_t marginal_start, uint64_t marginal_end, uint64_t end_overlap) {

    uint64_t alt_len = strlen(alt_token);
 
	struct lps substr;
	init_lps_offset(&substr, alt_token, alt_len, 0);
	lps_deepen(&substr, t_args->lcp_level);

    if (substr.size == 0) {
        // print first splitting node and link from reference graph	
        // print new node
        fprintf(t_args->out, "S\ts%d.%ld\t", t_args->thread_id, t_args->core_id_index);
        fwrite(chrom->seq+marginal_start, 1, start_loc-marginal_start, t_args->out);
        fprintf(t_args->out, "%s", alt_token);
        fwrite(chrom->seq+end_loc, 1, marginal_end-end_loc, t_args->out); 
        if (t_args->is_rgfa) {
            fprintf(t_args->out, "\tSN:Z:%d.%ld\tSO:i:%ld\tSR:i:1\n", t_args->thread_id, t_args->core_id_index, marginal_start); 
        } else {
            fprintf(t_args->out, "\n");
        }
        // print splitting link
        print_link(0, chrom->cores[latest_core_index].id, '+', t_args->thread_id, t_args->core_id_index, '+', 0, t_args->no_overlap, t_args->out);
        // print merging link
        print_link(t_args->thread_id, t_args->core_id_index, '+', 0, chrom->cores[first_core_after].id, '+', end_overlap, t_args->no_overlap, t_args->out);
        t_args->core_id_index++;
    } else {
        // print new node in between latest core before alternating core and first core in alternating core
        fprintf(t_args->out, "S\ts%d.%ld\t", t_args->thread_id, t_args->core_id_index);
        fwrite(chrom->seq+marginal_start, 1, start_loc-marginal_start, t_args->out);
        if (substr.cores[0].start) { // check if there is a string in alternating token that a cores in it don't cover, if do, merge it with a segment
            fwrite(alt_token, 1, substr.cores[0].start, t_args->out);
        }
        if (t_args->is_rgfa) {
            fprintf(t_args->out, "\tSN:Z:%d.%ld\tSO:i:%ld\tSR:i:1\n", t_args->thread_id, t_args->core_id_index, marginal_start); 
        } else {
            fprintf(t_args->out, "\n");
        }
        // print link between splitting segment with reference
        print_link(0, chrom->cores[latest_core_index].id, '+', t_args->thread_id, t_args->core_id_index, '+', 0, t_args->no_overlap, t_args->out);
        t_args->core_id_index++;

        // print first core as a segment
        char first_core_id[50];
        sprintf(first_core_id, "%d.%ld", t_args->thread_id, t_args->core_id_index);
        print_seq(t_args->thread_id, t_args->core_id_index, alt_token+substr.cores[0].start, substr.cores[0].end-substr.cores[0].start, first_core_id, substr.cores[0].start+start_loc, 0, t_args->is_rgfa, t_args->out);
        print_link(t_args->thread_id, t_args->core_id_index-1, '+', t_args->thread_id, t_args->core_id_index, '+', 0, t_args->no_overlap, t_args->out);
        t_args->core_id_index++;

        for (int i=1; i<substr.size; i++) {
            int start = substr.cores[i].start;
            int overlap = (int)(substr.cores[i-1].end-start);
            int offset = t_args->no_overlap ? overlap : 0;
            int core_len = substr.cores[i].end-start - offset;
            char temp_id[50];
            sprintf(temp_id, "%d.%ld", t_args->thread_id, t_args->core_id_index);
            print_seq(t_args->thread_id, t_args->core_id_index, alt_token+start+offset, core_len, temp_id, start_loc+start+offset, 1, t_args->is_rgfa, t_args->out);
            print_link(t_args->thread_id, t_args->core_id_index-1, '+', t_args->thread_id, t_args->core_id_index, '+', overlap, t_args->no_overlap, t_args->out);
            t_args->core_id_index++;
        }
        
        // create merging segment in between last core in alternating token and reference sequence.
        fprintf(t_args->out, "S\ts%d.%ld\t", t_args->thread_id, t_args->core_id_index);
        if (alt_len-substr.cores[substr.size-1].end) {
            fwrite(alt_token+substr.cores[substr.size-1].end, 1, alt_len-substr.cores[substr.size-1].end, t_args->out);
        }
        fwrite(chrom->seq+end_loc, 1, marginal_end-end_loc, t_args->out);
        if (t_args->is_rgfa) {
            char temp_id[50];
            sprintf(temp_id, "%d.%ld", t_args->thread_id, t_args->core_id_index);
            fprintf(t_args->out, "\tSN:Z:%s.%ld\tSO:i:%ld\tSR:i:1\n", temp_id, t_args->core_id_index, marginal_start); 
        } else {
            fprintf(t_args->out, "\n");
        }
        print_link(t_args->thread_id, t_args->core_id_index-1, '+', 0, chrom->cores[first_core_after].id, '+', end_overlap, t_args->no_overlap, t_args->out);
        t_args->core_id_index++;
    }

	free_lps(&substr);
}

void variate(struct t_arg *t_args, const struct chr *chrom, const char *org_seq, const char *alt_token, uint64_t start_loc) {

	// decide on boundaries
	uint64_t end_loc = start_loc+strlen(org_seq);
    uint64_t latest_core_index;
	uint64_t first_core_after;
	find_boundaries(start_loc, end_loc, chrom, &latest_core_index, &first_core_after);

	if (latest_core_index == 0 || chrom->cores[latest_core_index].end < chrom->cores[latest_core_index+1].start)  {
        t_args->failed_var_count += 1;
        pthread_mutex_lock(t_args->out_err_mutex);
        fprintf(t_args->out_err, "VARIATE-MARGIN-START:\tCHROM: %s,\tPOSITION: %ld,\tORG: %s,\tALT: %s,\tlatest_core_index: %ld\n", chrom->seq_name, start_loc, org_seq, alt_token, latest_core_index);
        fflush(t_args->out_err);
        pthread_mutex_unlock(t_args->out_err_mutex);
		return;
	}
	if (first_core_after+1 >= (uint64_t)chrom->cores_size || chrom->cores[first_core_after-1].end < chrom->cores[first_core_after].start) {
        t_args->failed_var_count += 1;
        pthread_mutex_lock(t_args->out_err_mutex);
        fprintf(t_args->out_err, "VARIATE-MARGIN-END:\tCHROM: %s,\tPOSITION: %ld,\tORG: %s,\tALT: %s,\tfirst_core_after: %ld,\tcores_size: %d\n", chrom->seq_name, start_loc, org_seq, alt_token, first_core_after, chrom->cores_size);
        fflush(t_args->out_err);
        pthread_mutex_unlock(t_args->out_err_mutex);
		return;
	}

    // get the chromosome
    uint64_t marginal_start = chrom->cores[latest_core_index].end; // start of variation
    uint64_t marginal_end = chrom->cores[first_core_after-1].end; // end of variation

    t_args->bubble_count += 1;

    // we need to make connection based on overlap, which is same size between first_core_after-1 and first_core_after
    uint64_t end_overlap = chrom->cores[first_core_after-1].end - chrom->cores[first_core_after].start;

    if (end_loc-start_loc<SV_LEN_BOUNDARY && strlen(alt_token)<SV_LEN_BOUNDARY) {
        variate_snp(t_args, chrom, alt_token, start_loc, end_loc, latest_core_index, first_core_after, marginal_start, marginal_end, end_overlap);
    } else {
        variate_sv(t_args, chrom, alt_token, start_loc, end_loc, latest_core_index, first_core_after, marginal_start, marginal_end, end_overlap);
    }

    return;
}