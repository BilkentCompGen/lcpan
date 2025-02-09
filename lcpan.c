#define STATS

/**
 * @file lcpan.cpp
 * @brief Implementation of a tool for constructing sequence graphs from genomic data.
 *
 * This file reads FASTA input and identifies its LCP (Locally Consistent Parsing) cores. 
 * It processes variations from a VCF file by locating the affected regions on the 
 * original sequence, identifying the LCP cores of those regions, aligning the original 
 * and varied cores, and creating bubbles to represent the differences. 
 * After applying all variations, it constructs a variation graph and outputs it in 
 * rGFA/GFA format, representing genomic sequences and their relationships for graph-based 
 * genome assembly and analysis.
 *
 */

#include "struct_def.h"
#include "opt_parser.h"
#include "fa_parser.h"
#include "vcf_parser.h"
#include "utils.h"

int main(int argc, char* argv[]) {

    struct opt_arg args;
    parse_opts(argc, argv, &args);

    LCP_INIT();

    struct ref_seq seqs; // sequence processed from fasta file
    
    read_fasta(&args, &seqs);

    FILE *out = fopen(args.gfa_path, "w");
    if (out == NULL) {
        fprintf(stderr, "Couldn't open output file %s\n", args.gfa_path);
        exit(EXIT_FAILURE);
    }

    print_ref_seq(&seqs, args.is_rgfa, args.no_overlap, out);

    FILE *out_err;
    if (args.prefix == NULL) {
        char out_err_filename[10];
        snprintf(out_err_filename, sizeof(out_err_filename), "lcpan.log");
        out_err = fopen(out_err_filename, "w");
    } else {
        char out_err_filename[strlen(args.prefix)+5];
        snprintf(out_err_filename, sizeof(out_err_filename), "%s.log", args.prefix);
        out_err = fopen(out_err_filename, "w");
    }
    
    if (out_err == NULL) {
        fprintf(stderr, "Couldn't open error log file\n");
        exit(EXIT_FAILURE);
    }
    fprintf(out_err, "%s\n", args.gfa_path);
    fprintf(out_err, "%d\n", args.thread_number);

    read_vcf(&args, &seqs, out_err);

    fclose(out);
    fclose(out_err);

    free_opt_arg(&args);
    free_ref_seq(&seqs);

    (void)(args.verbose && printf("[INFO] Total number of bubbles created: %d\n", args.bubble_count));
    (void)(args.verbose && printf("[INFO] Total number of invalid lines in the vcf file: %d\n", args.invalid_line_count));
    (void)(args.verbose && printf("[INFO] Total number of failed variations: %d\n", args.failed_var_count));
    
    return 0;
}