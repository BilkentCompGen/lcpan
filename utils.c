#include "utils.h"

int binary_search(uint64_t *arr, uint64_t size, uint64_t key) {
    int low = 0, high = size - 1;

    while (low <= high) {
        int mid = low + (high - low) / 2;
        if (arr[mid] == key) {
            return mid;
        }
        if (arr[mid] < key) {
            low = mid + 1;
        } else {
            high = mid - 1;
        }
    }
    return -1;
}

void quicksort(uint64_t *array, int low, int high) {
    if (low < high) {
        uint64_t pivot = array[high];
        int i = low - 1;

        for (int j = low; j < high; j++) {
            if (array[j] < pivot) {
                i++;
                uint64_t temp = array[i];
                array[i] = array[j];
                array[j] = temp;
            }
        }

        uint64_t temp = array[i+1];
        array[i+1] = array[high];
        array[high] = temp;

        quicksort(array, low, i);
        quicksort(array, i + 2, high);
    }
}

void print_seq3(uint64_t id, const char *seq1, int seq1_len, 
                             const char *seq2, int seq2_len,
                             const char *seq3, int seq3_len,
                             const char *seq_name, int start, int rank, int is_rgfa, FILE *out) {
	fprintf(out, "S\t%lu\t", id);
    fwrite(seq1, 1, seq1_len, out);
    fwrite(seq2, 1, seq2_len, out);
    fwrite(seq3, 1, seq3_len, out);
	if (is_rgfa) {
		fprintf(out, "\tSN:Z:%s\tSO:i:%d\tSR:i:%d\n", seq_name, start, rank); 
	} else {
		fprintf(out, "\n");
	}
}

void print_seq2(uint64_t id, const char *seq1, int seq1_len, 
                             const char *seq2, int seq2_len,
                             const char *seq_name, int start, int rank, int is_rgfa, FILE *out) {
	fprintf(out, "S\t%lu\t", id);
    fwrite(seq1, 1, seq1_len, out);
    fwrite(seq2, 1, seq2_len, out);
	if (is_rgfa) {
		fprintf(out, "\tSN:Z:%s\tSO:i:%d\tSR:i:%d\n", seq_name, start, rank); 
	} else {
		fprintf(out, "\n");
	}
}

void print_seq(uint64_t id, const char *seq, int seq_len, const char *seq_name, int start, int rank, int is_rgfa, FILE *out) {
	fprintf(out, "S\t%lu\t", id);
    fwrite(seq, 1, seq_len, out); 
	if (is_rgfa) {
		fprintf(out, "\tSN:Z:%s\tSO:i:%d\tSR:i:%d\n", seq_name, start, rank); 
	} else {
		fprintf(out, "\n");
	}
}

void print_seq3_vg(uint64_t id, const char *seq1, int seq1_len, 
                                const char *seq2, int seq2_len,
                                const char *seq3, int seq3_len,
                                const char *seq_name, int order, int start, int rank, int is_rgfa, FILE *out) {
	fprintf(out, "S\t%lu\t", id);
    fwrite(seq1, 1, seq1_len, out);
    fwrite(seq2, 1, seq2_len, out);
    fwrite(seq3, 1, seq3_len, out);
	if (is_rgfa) {
		fprintf(out, "\tSN:Z:%s.%d\tSO:i:%d\tSR:i:%d\n", seq_name, order, start, rank); 
	} else {
		fprintf(out, "\n");
	}
}

void print_seq2_vg(uint64_t id, const char *seq1, int seq1_len, 
                                const char *seq2, int seq2_len,
                                const char *seq_name, int order, int start, int rank, int is_rgfa, FILE *out) {
	fprintf(out, "S\t%lu\t", id);
    fwrite(seq1, 1, seq1_len, out);
    fwrite(seq2, 1, seq2_len, out);
	if (is_rgfa) {
		fprintf(out, "\tSN:Z:%s.%d\tSO:i:%d\tSR:i:%d\n", seq_name, order, start, rank); 
	} else {
		fprintf(out, "\n");
	}
}

void print_seq_vg(uint64_t id, const char *seq, int seq_len, const char *seq_name, int order, int start, int rank, int is_rgfa, FILE *out) {
	fprintf(out, "S\t%lu\t", id);
    fwrite(seq, 1, seq_len, out); 
	if (is_rgfa) {
		fprintf(out, "\tSN:Z:%s.%d\tSO:i:%d\tSR:i:%d\n", seq_name, order, start, rank); 
	} else {
		fprintf(out, "\n");
	}
}

void print_link(uint64_t id1, char sign1, uint64_t id2, char sign2, uint64_t overlap, FILE *out) {
    fprintf(out, "L\t%lu\t%c\t%lu\t%c\t%ldM\n", id1, sign1, id2, sign2, overlap);
}

void find_boundaries(uint64_t start_loc, uint64_t end_loc, const struct chr *chrom, uint64_t start_index, uint64_t *latest_core_index, uint64_t *first_core_after) {

    const struct simple_core *cores = chrom->cores;
    uint64_t cores_size = (uint64_t)chrom->cores_size;
    uint64_t core_idx_1 = cores[start_index].end <= start_loc ? start_index : 0;
    int left, mid, right;

    left = core_idx_1;
    right = core_idx_1+5 < cores_size && end_loc < cores[core_idx_1+5].start ? core_idx_1+5 : cores_size;
    while (left < right) {
        mid = left + (right - left) / 2;
        if (cores[mid].end <= start_loc) {
            core_idx_1 = mid;
            left = mid + 1;
        } else {
            right = mid;
        }
    }

    while (core_idx_1 && start_loc < cores[core_idx_1].end ) {
        core_idx_1--;
    }

    *latest_core_index = core_idx_1;
    *first_core_after = cores_size;
	for (uint64_t j=core_idx_1; j<cores_size; j++) {
        if (end_loc <= cores[j].start) {
            *first_core_after = j;
			return;
        }
    }
}

void refine_seqs(struct ref_seq *seqs, int no_overlap) {
    if (no_overlap) {
        for (int i=0; i<seqs->size; i++) {          
            for (int j=1; j<seqs->chrs[i].cores_size; j++) {
                if (seqs->chrs[i].cores[j].start < seqs->chrs[i].cores[j-1].end) {
                    seqs->chrs[i].cores[j].start = seqs->chrs[i].cores[j-1].end;
                }
            }
        }
    }
}

void refine_seq(struct lps *str, int no_overlap) {
    if (no_overlap) {
        for (int j=1; j<str->size; j++) {
            if (str->cores[j].start < str->cores[j-1].end) {
                str->cores[j].start = str->cores[j-1].end;
            }
        }
    }
}

void print_path(const struct ref_seq *seqs, FILE *out) {
    for (int i=0; i<seqs->size; i++) {
        if (seqs->chrs[i].cores_size) {
            const struct chr *chrom = seqs->chrs + i;

            // print Path (P)
            fprintf(out, "P\t%s\t", chrom->seq_name);

            // print first core
            if (chrom->ids != NULL && chrom->ids[0] != NULL) { // if first one is not NULL
                fprintf(out, "%lu+", chrom->ids[0][0]);
                int index = 1;
                while (chrom->ids[0][index]) {
                    fprintf(out, ",%lu+", chrom->ids[0][index++]);
                }
                fprintf(out, ",%lu+", chrom->cores[0].id);
            } else {
                fprintf(out, "%lu+", chrom->cores[0].id);
            }
            
            // print rest
            if (chrom->ids != NULL) {
                for (int j=1; j<chrom->cores_size; j++) {
                    if (chrom->ids[j] != NULL) {
                        int index = 0;
                        while (chrom->ids[j][index]) {
                            fprintf(out, ",%lu+", chrom->ids[j][index++]);
                        }
                    }
                    fprintf(out, ",%lu+", chrom->cores[j].id);
                }
            } else {
                for (int j=1; j<chrom->cores_size; j++) {
                    fprintf(out, ",%lu+", chrom->cores[j].id);
                }
            }

            // print cigar
            fprintf(out, "\t*\n");
        }
    }
}
