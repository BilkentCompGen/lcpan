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

void process_chrom_vg(char *sequence, uint64_t seq_size, int lcp_level, struct chr *chrom, uint64_t *core_id_index) {
    uint64_t id = *core_id_index;

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
            chrom->cores[last_core_index].id = id;
            chrom->cores[last_core_index].start = str.cores[i].start;
            chrom->cores[last_core_index].end = str.cores[i].end;
            id++;
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

    *core_id_index = id;
}

void process_chrom_ldbg(char *sequence, uint64_t seq_size, int lcp_level, struct chr *chrom) {
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

    printf("[INFO] Started processing reference file.\n");

    // read index (fai) file to get chromosome count and allocate space for chromosomes for later processing
    char line[1024];
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

    while (fgets(line, sizeof(line), idx) != NULL) {
        char *name, *length;
        
        // assign name
        name = strtok(line, "\t");
        uint64_t name_len = strlen(name);
        seqs->chrs[chrom_index].seq_name = (char *)malloc(name_len);
        memcpy(seqs->chrs[chrom_index].seq_name, name, name_len);
        
        // assign size and allocate in memory
        length = strtok(NULL, "\t");
        seqs->chrs[chrom_index].seq_size = strtol(length, NULL, 10);

        seqs->chrs[chrom_index].seq = (char *)malloc(seqs->chrs[chrom_index].seq_size);
        if (seqs->chrs[chrom_index].seq == NULL) {
            fprintf(stderr, "REF: Couldn't allocate memory to chromosome string.\n");
            exit(EXIT_FAILURE);
        }
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
    while (fgets(line, sizeof(line), ref)) {

        line[strcspn(line, "\n")] = '\0';

        if (line[0] == '>') {
            if (sequence_size != 0) {
                if (args->program == VG) {
                    process_chrom_vg(seqs->chrs[index].seq, sequence_size, args->lcp_level, &(seqs->chrs[index]), &(args->core_id_index));
                } else if (args->program == LDBG) {
                    process_chrom_ldbg(seqs->chrs[index].seq, sequence_size, args->lcp_level, &(seqs->chrs[index]));
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
        if (args->program == VG) {
            process_chrom_vg(seqs->chrs[index].seq, sequence_size, args->lcp_level, &(seqs->chrs[index]), &(args->core_id_index));
        } else if (args->program == LDBG) {
            process_chrom_ldbg(seqs->chrs[index].seq, sequence_size, args->lcp_level, &(seqs->chrs[index]));
        }
        index++;
    }

    fclose(ref);

    printf("[INFO] Input file processing ended.\n");
}
