#include "lbdg.h"

void lbdg_print_ref_seq(struct ref_seq *seqs, FILE *out) {

    printf("[INFO] Processing reference...\n");

    fprintf(out, "H\tVN:Z:1.1\n");

    uint64_t total_cores = 0;
    for (int i=0; i<seqs->size; i++) total_cores += seqs->chrs[i].cores_size;

    uint64_t *segment_ids = (uint64_t *)malloc(total_cores * sizeof(uint64_t));
    uint64_t index = 0;

    for (int i=0; i<seqs->size; i++) {
        const struct chr chrom = seqs->chrs[i];
        for (int j=0; j<chrom.cores_size; j++) segment_ids[index++] = chrom.cores[j].id;
    }

    quicksort(segment_ids, 0, index-1);

    uint64_t size=0, iter=1;

    while (iter<index) {
        if (segment_ids[size] != segment_ids[iter]) {
            size++;
            segment_ids[size] = segment_ids[iter];
        }
        iter++;
    }

    size++; // total distinct count

    uint64_t *bit_arr = (uint64_t *)calloc(size, sizeof(uint64_t));

    for (int i=0; i<seqs->size; i++) {
		struct simple_core *cores = seqs->chrs[i].cores;
		const char *seq = seqs->chrs[i].seq;

        if (seqs->chrs[i].cores_size) {
            int index = binary_search(segment_ids, size, cores[0].id);
            if (index != -1) {
                uint64_t bit = (1ULL << (index % 64));
                if (!(bit_arr[index / 64] & bit)) {
                    const struct simple_core *curr_core = &(seqs->chrs[i].cores[0]);
                    uint64_t curr_start = curr_core->start;
                    uint64_t curr_end = curr_core->end;
                    int seq_len = (int)(curr_end-curr_start);

                    fprintf(out, "S\t%lu\t", curr_core->id);
                    fwrite(seq+curr_start, 1, seq_len, out); 
                    fprintf(out, "\n");

                    bit_arr[index / 64] |= bit;
                }
            }
        }
		
		for (int j=1; j<seqs->chrs[i].cores_size; j++) {
            // binary search the index of curr_core->id in segmet_ids and if in that index
            int index = binary_search(segment_ids, size, cores[j].id);
            if (index != -1) {
                const struct simple_core *curr_core = &(seqs->chrs[i].cores[j]);
                const struct simple_core *prev_core = &(seqs->chrs[i].cores[j-1]);
                uint64_t curr_start = curr_core->start;

                uint64_t bit = (1ULL << (index % 64));
                if (!(bit_arr[index / 64] & bit)) {
                    
                    uint64_t curr_end = curr_core->end;
                    int seq_len = (int)(curr_end-curr_start);

                    fprintf(out, "S\t%lu\t", curr_core->id);
                    fwrite(seq+curr_start, 1, seq_len, out); 
                    fprintf(out, "\n");

                    bit_arr[index / 64] |= bit;
                }

                int overlap = (int)(prev_core->end-curr_start);
                if (0 < overlap) {
                    fprintf(out, "L\t%lu\t+\t%lu\t+\t%dM\n", prev_core->id, curr_core->id, overlap);
                }
            }
		}
	}

    free(bit_arr);
}