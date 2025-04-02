#include "fa_parser.h"

uint64_t MurmurHash3_32(const void *key, int len) {
    const uint8_t *data = (const uint8_t *)key;
    const int nblocks = len / 4;

    uint32_t h1 = 42;

    const uint32_t c1 = 0xcc9e2d51;
    const uint32_t c2 = 0x1b873593;

    // Body: Process blocks of 4 bytes at a time
    const uint32_t *blocks = (const uint32_t *)(data + nblocks * 4);

    for (int i = -nblocks; i; i++) {
        uint32_t k1 = blocks[i];

        k1 *= c1;
        k1 = (k1 << 15) | (k1 >> (32 - 15));
        k1 *= c2;

        h1 ^= k1;
        h1 = (h1 << 15) | (h1 >> (32 - 15));
        h1 = h1 * 5 + 0xe6546b64;
    }

    // Tail: Process remaining bytes
    const uint8_t *tail = (const uint8_t *)(data + nblocks * 4);

    uint32_t k1 = 0;

    switch (len & 3) {
    case 3:
        k1 ^= tail[2] << 16;
        break;
    case 2:
        k1 ^= tail[1] << 8;
        break;
    case 1:
        k1 ^= tail[0];
        k1 *= c1;
        k1 = (k1 << 15) | (k1 >> (32 - 15));
        k1 *= c2;
        h1 ^= k1;
    }

    // Finalization: Mix the hash to ensure the last few bits are fully mixed
    h1 ^= len;

    /* fmix32 */
    h1 ^= h1 >> 16;
    h1 *= 0x85ebca6b;
    h1 ^= h1 >> 13;
    h1 *= 0xc2b2ae35;
    h1 ^= h1 >> 16;
    return (uint64_t)h1;
}

void free_ref_seq(struct ref_seq *seqs) {
	if (seqs->size) {
		for (int i=0; i<seqs->size; i++) {
			free(seqs->chrs[i].seq_name);
			free(seqs->chrs[i].seq);
            if (seqs->chrs[i].cores_size) {
			    free(seqs->chrs[i].cores);
            }
            
            if (seqs->chrs[i].ids) {
                for (int j=0; j<seqs->chrs[i].cores_size; j++) {
                    if (seqs->chrs[i].ids[j] != NULL) free(seqs->chrs[i].ids[j]);
                }
                free(seqs->chrs[i].ids);
            }
		}
		free(seqs->chrs);
		seqs->size = 0;
	}
}

void vgx_process_chrom(char *sequence, uint64_t seq_size, int lcp_level, int skip_masked, struct chr *chrom, uint64_t *core_id_index) {
    uint64_t id = *core_id_index;

    uint64_t estimated_core_size = (uint64_t)(seq_size / pow(1.5, lcp_level));
    uint64_t estimated_core_length = (uint64_t)(3 * pow(2, lcp_level-1));
    chrom->cores_size = 0;

    if (estimated_core_size == 0) {
        chrom->cores = NULL;
        return;
    }
    
    chrom->cores = (struct simple_core*)malloc(estimated_core_size * sizeof(struct simple_core));

    uint64_t index = 0;
    uint64_t last_core_index = 0;

    while (index < seq_size) {
        {
            uint64_t temp_index = index;
            
            while (index < seq_size && sequence[index] == 'N') {
                index++;
            }

            if (!skip_masked) {
                if (temp_index != index) {
                    uint64_t i = temp_index;
                    while (i+estimated_core_length < index) {
                        chrom->cores[last_core_index].id = id;
                        chrom->cores[last_core_index].start = i;
                        chrom->cores[last_core_index].end = i+estimated_core_length;
                        id++;
                        last_core_index++;
                        i += estimated_core_length;
                    }
                    chrom->cores[last_core_index].id = id;
                    chrom->cores[last_core_index].start = i;
                    chrom->cores[last_core_index].end = index;
                    id++;
                    last_core_index++;
                }
            }
        }

        uint64_t end = index;
        
        while (end < seq_size && sequence[end] != 'N') {
            end++;
        }

        struct lps str;
        init_lps_offset(&str, sequence+index, end-index, index);
        lps_deepen(&str, lcp_level);

        if (str.size) {
            if (str.cores[0].start != index) {
                chrom->cores[last_core_index].id = id;
                chrom->cores[last_core_index].start = index;
                chrom->cores[last_core_index].end = str.cores[0].start;
                id++;
                last_core_index++;                
            }

            for (int i=0; i<str.size; i++) {
                chrom->cores[last_core_index].id = id;
                chrom->cores[last_core_index].start = str.cores[i].start;
                chrom->cores[last_core_index].end = str.cores[i].end;
                id++;
                last_core_index++;
            }
    
            if (str.cores[str.size-1].end != end) {
                chrom->cores[last_core_index].id = id;
                chrom->cores[last_core_index].start = str.cores[str.size-1].end;
                chrom->cores[last_core_index].end = end;
                id++;
                last_core_index++;                
            }
            
            chrom->cores_size = last_core_index;
        }

        index = end;
        free_lps(&str); 
    }

    if (chrom->cores_size) {
        struct simple_core *temp = (struct simple_core*)realloc(chrom->cores, chrom->cores_size * sizeof(struct simple_core));
        if (temp != NULL) {
            chrom->cores = temp;
        }
    } else {
        free(chrom->cores);
        chrom->cores = NULL;
    }

    *core_id_index = id;
}

void lbdg_process_chrom(char *sequence, uint64_t seq_size, int lcp_level, struct chr *chrom) {
    uint64_t estimated_core_size = (int)(seq_size / pow(1.5, lcp_level));
    chrom->cores_size = 0;

    if (estimated_core_size == 0) {
        chrom->cores = NULL;
        return;
    }
    
    chrom->cores = (struct simple_core*)malloc(estimated_core_size * sizeof(struct simple_core));

    uint64_t index = 0;
    uint64_t last_core_index = 0;

    while (index < seq_size) {
        while (index < seq_size && sequence[index] == 'N') {
            index++;
        }
        uint64_t end = index;
        
        while (end < seq_size && sequence[end] != 'N') {
            end++;
        }

        struct lps str;
        init_lps_offset(&str, sequence+index, end-index, index);
        lps_deepen(&str, lcp_level);

        for (int i=0; i<str.size; i++) {
            chrom->cores[last_core_index].id = (MurmurHash3_32(sequence + str.cores[i].start, (int)(str.cores[i].end-str.cores[i].start)) << 32) | str.cores[i].label;
            chrom->cores[last_core_index].start = str.cores[i].start;
            chrom->cores[last_core_index].end = str.cores[i].end;
            last_core_index++;
        }

        chrom->cores_size = last_core_index;
        index = end;

        free_lps(&str); 
    }

    struct simple_core *temp = (struct simple_core*)realloc(chrom->cores, chrom->cores_size * sizeof(struct simple_core));

    if (temp != NULL) {
        chrom->cores = temp;
    }
}

void read_fasta(struct opt_arg *args, struct ref_seq *seqs) {

    printf("[INFO] Processing reference...\n");

    // read index (fai) file to get chromosome count and allocate space for chromosomes for later processing
    int line_size = 1024;
    char line[line_size];
    FILE *idx = fopen(args->fasta_fai_path, "r");
    if (idx == NULL) {
        fprintf(stderr, "REF: Couldn't open file %s\n", args->fasta_fai_path);
        exit(EXIT_FAILURE);
    }
    
    int chrom_index = 0;

    while (fgets(line, sizeof(line), idx) != NULL) {
        chrom_index++;
    }

    if (chrom_index == 0) {
        fprintf(stderr, "REF: Index file is empty.\n");
        exit(EXIT_FAILURE);
    }

    seqs->size = chrom_index;
    seqs->chrs = (struct chr *)malloc(chrom_index*sizeof(struct chr));
    if (seqs->chrs == NULL) {
        fprintf(stderr, "REF: Couldn't allocate memory to ref sequences\n");
        exit(EXIT_FAILURE);
    }

    rewind(idx);
    chrom_index = 0;
    uint64_t global_index = 0;

    while (fgets(line, sizeof(line), idx) != NULL) {
        char *name, *length;
        
        // assign name
        char *saveptr;
        name = strtok_r(line, "\t", &saveptr);
        seqs->chrs[chrom_index].global_index = global_index;
        uint64_t name_len = strlen(name);
        seqs->chrs[chrom_index].seq_name = (char *)malloc(name_len+1);
        memcpy(seqs->chrs[chrom_index].seq_name, name, name_len);
        seqs->chrs[chrom_index].seq_name[name_len] = '\0';

        // assign size and allocate in memory
        length = strtok_r(NULL, "\t", &saveptr);
        seqs->chrs[chrom_index].seq_size = strtol(length, NULL, 10);
        global_index += seqs->chrs[chrom_index].seq_size;

        seqs->chrs[chrom_index].seq = (char *)malloc(seqs->chrs[chrom_index].seq_size+1);
        if (seqs->chrs[chrom_index].seq == NULL) {
            fprintf(stderr, "REF: Couldn't allocate memory to chromosome string.\n");
            exit(EXIT_FAILURE);
        }
        seqs->chrs[chrom_index].seq[seqs->chrs[chrom_index].seq_size] = '\0';
        seqs->chrs[chrom_index].ids = NULL;
        chrom_index++;
    }

    fclose(idx);

    // move to read reference file
    FILE *ref = fopen(args->fasta_path, "r");
    if (ref == NULL) {
        fprintf(stderr, "REF: Couldn't open file %s\n", args->fasta_path);
        exit(EXIT_FAILURE);
    }

    // make necesarry declaretions and initialization
    uint64_t sequence_size = 0;
    int index = 0;

    // process reference file
    while (fgets(line, line_size, ref)) {

        line[strcspn(line, "\n")] = '\0';

        if (line[0] == '>') {
            if (sequence_size != 0) {
                if (args->program == VG || args->program == VGX) {
                    vgx_process_chrom(seqs->chrs[index].seq, sequence_size, args->lcp_level, args->skip_masked, &(seqs->chrs[index]), &(args->core_id_index));
                } else if (args->program == LBDG) {
                    lbdg_process_chrom(seqs->chrs[index].seq, sequence_size, args->lcp_level, &(seqs->chrs[index]));
                }
                sequence_size = 0;
                index++;
            }
        } else {
            uint64_t line_len = strlen(line);
            memcpy(seqs->chrs[index].seq + sequence_size, line, line_len);
            sequence_size += line_len;
        }
    }

    if (sequence_size != 0) {
        if (args->program == VG || args->program == VGX) {
            vgx_process_chrom(seqs->chrs[index].seq, sequence_size, args->lcp_level, args->skip_masked, &(seqs->chrs[index]), &(args->core_id_index));
        } else if (args->program == LBDG) {
            lbdg_process_chrom(seqs->chrs[index].seq, sequence_size, args->lcp_level, &(seqs->chrs[index]));
        }
        index++;
    }

    fclose(ref);
}

void print_ref_seqs(const struct ref_seq *seqs, int is_rgfa, FILE *out) {

    printf("[INFO] Printing reference...\n");

    fprintf(out, "H\tVN:Z:1.1\n");

	// iterate through each chromosome
	for (int i=0; i<seqs->size; i++) {
		if (seqs->chrs[i].cores_size) {
            const char *seq_name = seqs->chrs[i].seq_name;
            const char *seq = seqs->chrs[i].seq;
            
            {
                const struct simple_core *curr_core = &(seqs->chrs[i].cores[0]);
                uint64_t start = curr_core->start;
                print_seq(curr_core->id, seq+start, (int)(curr_core->end - start), seq_name, start, 0, is_rgfa, out);
            }
            
            for (int j=1; j<seqs->chrs[i].cores_size; j++) {
                const struct simple_core *curr_core = &(seqs->chrs[i].cores[j]);
                const struct simple_core *prev_core = &(seqs->chrs[i].cores[j-1]);
                uint64_t curr_start = curr_core->start;
                int seq_len = (int)(curr_core->end - curr_start);
                int overlap = (int)(prev_core->end - curr_start);
                
                // there might be graps ('N') in genome, hence
                if (overlap < 0) {
                    overlap = 0;
                }

                print_seq(curr_core->id, seq+curr_start, seq_len, seq_name, curr_start, 0, is_rgfa, out);
                print_link(prev_core->id, '+', curr_core->id, '+', overlap, out);
            }

            // Print Path (P)
            fprintf(out, "P\t%s\t", seq_name);
            fprintf(out, "%lu+", seqs->chrs[i].cores[0].id);
            for (int j=1; j<seqs->chrs[i].cores_size; j++) {
                fprintf(out, ",%lu+", seqs->chrs[i].cores[j].id);
            }
            fprintf(out, "\t*\n");
        }
	}
}