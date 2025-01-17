#include "vcf_parser.h"

long count = 0;

void read_vcf(struct opt_arg *args, struct ref_seq *seqs, int* failed_var_count, int* invalid_line_count, FILE *out, FILE *out_err) {
    char line[1024];
    
    FILE *file = fopen(args->vcf_path, "r");
    if (file == NULL) {
        fprintf(stderr, "Couldn't open file %s\n", args->vcf_path);
        exit(EXIT_FAILURE);
    }

    long last_position = 0;

    while (fgets(line, sizeof(line), file) != NULL) {
        line[strlen(line)-1] = '\0';
        // char *tag, *value;

        
        
        if (line[0] == '#' && line[1] == '#') {
            // tag = strtok(line, "="); 
            // if (tag == NULL) {
            //     continue; // Skip if no '=' found
            // }
            // tag += 2; // Remove "##" 
            // value = strtok(NULL, "\n"); 

            // if (strncmp(tag, "INFO", 4) == 0) {
            //     // handle ##INFO=<line> 
            //     printf("INFO: %s\n", value);
            //     value++;                        // remove <
            //     value[strlen(value)-1] = '\0';  // remove >

            //     // parse INFO tags
            //     char *info_tag, *info_value;
            //     info_tag = strtok(value, "=");
            //     while (info_tag != NULL) {
            //         info_value = strtok(NULL, ",");
            //         printf("tag: %s, value: %s\n", info_tag, info_value);
            //         info_tag = strtok(NULL, "=");
            //     }
            // } else if (strncmp(tag, "FILTER", 6) == 0) {
            //     // // handle ##FILTER=<line>
            //     printf("FILTER: %s\n", value);
            //     value++;                        // remove <
            //     value[strlen(value)-1] = '\0';  // remove >

            //     // parse FILTER tags
            //     char *id, *description;
            //     id = strtok(value, "=");
            //     id = strtok(NULL, ",");
            //     description = strtok(NULL, "=");
            //     description = strtok(NULL, ",");
            //     description++;
            //     description[strlen(description)-1] = '\0';
            //     printf("ID: %s, Description: %s\n", id, description);
            // } else if (strncmp(tag, "FORMAT", 6) == 0) {
            //     // Handle ##FORMAT=<line>
            //     printf("FORMAT: %s\n", value); 
            // }
        } else if (line[0] == '#') { // header line
            // printf("HEADER: %s\n", line);
            // char *format_ptr = strstr(line, "FORMAT");
            // if (format_ptr == NULL) {
            //     fprintf(stderr, "Error: 'FORMAT' field not found in header.\n");
            //     return;
            // }
        } else {
            break;
        }

        last_position = ftell(file);
    }

    fseek(file, last_position, SEEK_SET);

    while (fgets(line, sizeof(line), file) != NULL) {
        line[strlen(line)-1] = '\0';
        // printf("%s\n", line);
        
        char *chrom, *index, *seq, *alt;
        long offset;
    
        // split the line by tab characters
        if (line == NULL) {
            *invalid_line_count = (*invalid_line_count) + 1;
            continue;
        }
        chrom = strtok(line, "\t"); // get chromosome name
        
        index = strtok(NULL, "\t"); // get index

        if (index == NULL) {
            *invalid_line_count = (*invalid_line_count) + 1;
            continue;
        }
        offset = strtol(index, NULL, 10);
        strtok(NULL, "\t");         // skip ID
        seq = strtok(NULL, "\t");   // get sequence
        alt = strtok(NULL, "\t");   // get ALT alleles

        // printf("chrom: %s, index %ld, seq: %s\n", chrom, offset, seq);

        int chrom_index = -1;
        for (int i=0; i<seqs->size; i++) {
            if (strcmp(chrom, seqs->chrs[i].seq_name) == 0) {
                chrom_index = i;
                break;
            }
        }

        if (chrom_index == -1) {
            fprintf(stderr, "Couldn't locate chrom %s from VCF in referefence\n", chrom);
            *invalid_line_count = (*invalid_line_count) + 1;
            continue;
        }

        char *alt_token = strtok(alt, ","); // split ALT alleles by comma
        while (alt_token != NULL) {
            // printf("- %s - len: %ld\n", alt_token, strlen(alt_token));

            count++;
            printf("#%ld, chrom: %s, index %ld, seq: %s, alt: %s\n", count, chrom, offset, seq, alt_token);

            // handle variation            
            variate(&(seqs->chrs[chrom_index]), seq, alt_token, offset, args->lcp_level, &(args->core_id_index), failed_var_count, out, out_err);

            alt_token = strtok(NULL, ",");
        }
    }

    fclose(file);
}