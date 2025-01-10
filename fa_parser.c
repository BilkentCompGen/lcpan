#include "fa_parser.h"

void process_chrom(char *sequence, size_t seq_size, int lcp_level, struct chr *chrom, uint64_t *core_id_index) {
    uint64_t id = *core_id_index;

    struct lps str;
    init_lps(&str, sequence, seq_size);
    lps_deepen(&str, lcp_level);

    chrom->cores_size = str.size;
    chrom->cores = (struct simple_core *)malloc((chrom->cores_size)*sizeof(struct simple_core));

    for (int i=0; i<str.size; i++) {
        chrom->cores[i].id = id;
        chrom->cores[i].start = str.cores[i].start;
        chrom->cores[i].end = str.cores[i].end;
        chrom->cores[i].label = str.cores[i].label;
        id++;
    }

    *core_id_index = id;
    free_lps(&str);
}

void read_fasta(struct opt_arg *args, struct ref_seq *seqs) {

    // read index (fai) file to get chromosome count and allocate space for chromosomes for later processing
    char line[1024];
    FILE *idx = fopen(args->fasta_fai_path, "r");
    if (idx == NULL) {
        fprintf(stderr, "Couldn't open file %s\n", args->fasta_fai_path);
        exit(EXIT_FAILURE);
    }
    
    int chrom_index = 0;

    while (fgets(line, sizeof(line), idx) != NULL) {
        chrom_index++;
    }

    if (chrom_index == 0) {
        fprintf(stderr, "Index file is empty.\n");
        exit(EXIT_FAILURE);
    }

    seqs->size = chrom_index;
    seqs->chrs = (struct chr *)malloc(chrom_index*sizeof(struct chr));
    if (seqs->chrs == NULL) {
        fprintf(stderr, "Couldn't allocate memory to ref sequences\n");
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
            fprintf(stderr, "Couldn't allocate memory to chromosome string.\n");
            exit(EXIT_FAILURE);
        }
        chrom_index++;
    }

    fclose(idx);

    // move to read reference file
    FILE *ref = fopen(args->fasta_path, "r");
    if (ref == NULL) {
        fprintf(stderr, "Couldn't open file %s\n", args->fasta_path);
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
                process_chrom(seqs->chrs[index].seq, sequence_size, args->lcp_level, &(seqs->chrs[index]), &(args->core_id_index));
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
        process_chrom(seqs->chrs[index].seq, sequence_size, args->lcp_level, &(seqs->chrs[index]), &(args->core_id_index));
        index++;
    }

    fclose(ref);
}
