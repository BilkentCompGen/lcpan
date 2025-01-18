#define STATS

/**
 * @file vg_lcp.cpp
 * @brief Implementation of a tool for constructing sequence graphs from genomic data.
 *
 * This file reads FASTA input and identifies its LCP (Locally Consistent Parsing) cores. 
 * It processes variations from a VCF file by locating the affected regions on the 
 * original sequence, identifying the LCP cores of those regions, aligning the original 
 * and varied cores, and creating bubbles to represent the differences. 
 * After applying all variations, it constructs a variation graph and outputs it in 
 * rGFA format, representing genomic sequences and their relationships for graph-based 
 * genome assembly and analysis.
 *
 * Structures used:
 * - `segment`: Represents a segment in the rGFA format.
 * - `g_link`: Represents a graph link in the rGFA format.
 * - `core_node`: Represents a node in the sequence graph with pointers to other nodes 
 *   and its LCP core value.
 *
 * Utility functions:
 * - `core_to_seq()`: Converts a core string representation to a nucleotide sequence.
 * - `dfs()`: Performs a depth-first traversal of the sequence graph to generate graph links.
 * - `find_boundaries()`: Determines the boundaries of a variation in the sequence graph.
 * - `align_variation()`: Aligns a variation against the sequence graph to identify differences.
 * - `construct_variated_seq()`: Constructs a modified sequence incorporating a variation.
 * - `variate()`: Applies variations to the sequence graph, creating bubble structures.
 *
 * Command-line options:
 * - `-f`: Specifies the input FASTA file path.
 * - `-v`: Specifies the input VCF file path.
 * - `-r`: Specifies the output rGFA file path.
 * - `-l`: Specifies the LCP level for parsing.
 */

#include "struct_def.h"
#include "opt_parser.h"
#include "fa_parser.h"
#include "vcf_parser.h"
#include "utils.h"

/**
 * @brief Main function to process the input, construct the sequence graph, and apply variations.
 * @param argc The number of command-line arguments.
 * @param argv The array of command-line arguments.
 * @return Exit status of the program (0 for success, non-zero for errors).
 */
int main(int argc, char* argv[]) {

    struct opt_arg args;
    parse_opts(argc, argv, &args);

    LCP_INIT();

    struct ref_seq seqs; // sequence processed from fasta file
    
    read_fasta(&args, &seqs);

    FILE *out = fopen(args.rgfa_path, "w");
    if (out == NULL) {
        fprintf(stderr, "Couldn't open output file %s\n", args.rgfa_path);
        exit(EXIT_FAILURE);
    }

    print_ref_seq(&seqs, out);

    FILE *out_err = fopen("out.err", "w");
    if (out_err == NULL) {
        fprintf(stderr, "Couldn't open output file %s\n", "out_err");
        exit(EXIT_FAILURE);
    }

    int failed_var_count = 0;
    int invalid_line_count = 0;
			 

    read_vcf(&args, &seqs, &failed_var_count, &invalid_line_count, out, out_err);

    fclose(out);
    fclose(out_err);

    free_opt_arg(&args);
    free_ref_seq(&seqs);

    printf("Total number of failed variations: %d\n", failed_var_count);
    printf("Total number of invalid lines in the vcf file: %d\n", invalid_line_count);
    
    return 0;
}