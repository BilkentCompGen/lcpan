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
#include "utils.h"
#include "opt_parser.h"
#include "fa_parser.h"
#include "vg.h"
#include "vgx.h"
#include "lbdg.h"

int main(int argc, char* argv[]) {

    struct opt_arg args;
    parse_opts(argc, argv, &args);

    LCP_INIT();

    struct ref_seq seqs; // sequence processed from fasta file

    read_fasta(&args, &seqs);

    FILE *gfa_out;

    switch (args.program) {
    
    case VG:
        refine_seqs(&seqs, args.no_overlap);
        vg_read_vcf(&args, &seqs);
        break;
    case VGX:
        refine_seqs(&seqs, args.no_overlap);
        gfa_out = fopen(args.gfa_path, "w");
        if (gfa_out == NULL) {
            fprintf(stderr, "Couldn't open output file %s\n", args.gfa_path);
            exit(EXIT_FAILURE);
        }
        print_ref_seqs(&seqs, args.is_rgfa, gfa_out);
        vgx_read_vcf(&args, &seqs);
        (void)(args.verbose && printf("[INFO] Total number of bubbles created: %d\n", args.bubble_count));
        (void)(args.verbose && printf("[INFO] Total number of invalid lines in the vcf file: %d\n", args.invalid_line_count));
        (void)(args.verbose && printf("[INFO] Total number of failed variations: %d\n", args.failed_var_count));
        fclose(gfa_out);
        break;
    case LBDG:
        gfa_out = fopen(args.gfa_path, "w");
        if (gfa_out == NULL) {
            fprintf(stderr, "Couldn't open output file %s\n", args.gfa_path);
            exit(EXIT_FAILURE);
        }
        lbdg_print_ref_seq(&seqs, gfa_out);
        fclose(gfa_out);
        break;
    default:
        fprintf(stderr, "Invalid program mode provided.\n");
    }
    
    free_opt_arg(&args);
    free_ref_seq(&seqs);

    return 0;
}