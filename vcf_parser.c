#include "vcf_parser.h"

void read_vcf(struct opt_arg *args, struct ref_seq *seqs, int* failed_var_count, int* invalid_line_count, int* bubble_count, FILE *out, FILE *out_err) {
    
    printf("[INFO] Started processing variation file.\n");
    
    uint64_t current_size = 1048576;
    char *line = (char *)malloc(current_size);
    
    FILE *file = fopen(args->vcf_path, "r");
    if (file == NULL) {
        fprintf(out_err, "VCF: Couldn't open file %s\n", args->vcf_path);
        exit(EXIT_FAILURE);
    }

    while (fgets(line, current_size, file) != NULL) {
        while (strlen(line) == current_size - 1 && line[current_size - 2] != '\n') {
            long current_position = ftell(file);
            fseek(file, current_position - (current_size - 1), SEEK_SET); // last char is '\0'
            
            current_size *= 2;
            line = (char *)realloc(line, current_size);
            continue; // retry to read line
        }

        line[strlen(line)-1] = '\0';
        
        if (line[0] == '#')
            continue;
        
        char *chrom, *index, *seq, *alt;
        long offset;

        // split the line by tab characters
        if (line == NULL) {
            *invalid_line_count = (*invalid_line_count) + 1;
            fprintf(out_err, "VCF: line: %s\n", line);
            continue;
        }

        // get chromosome name
        chrom = strtok(line, "\t"); 

        // get index
        index = strtok(NULL, "\t"); 
        if (index == NULL) {
            *invalid_line_count = (*invalid_line_count) + 1;
            fprintf(out_err, "VCF: no index at line: %s\n", line);
            continue;
        }

        // get offset
        offset = strtol(index, NULL, 10);
        strtok(NULL, "\t");         // skip ID

        // get sequence
        seq = strtok(NULL, "\t"); 

        // get ALT alleles
        alt = strtok(NULL, "\t");

        int chrom_index = -1;
        for (int i=0; i<seqs->size; i++) {
            if (strcmp(chrom, seqs->chrs[i].seq_name) == 0) {
                chrom_index = i;
                break;
            }
        }

        if (chrom_index == -1) {
            fprintf(out_err, "VCF: Couldn't locate chrom %s from VCF in reference\n", chrom);
            *invalid_line_count = (*invalid_line_count) + 1;
            continue;
        }

        char *alt_token = strtok(alt, ","); // split ALT alleles by comma
        
        while (alt_token != NULL) {
            variate(&(seqs->chrs[chrom_index]), seq, alt_token, offset, args->lcp_level, &(args->core_id_index), failed_var_count, bubble_count, args->is_rgfa, out, out_err);
            alt_token = strtok(NULL, ",");
        }
    }

    fclose(file);

    printf("[INFO] Ended processing variation file.\n");
}