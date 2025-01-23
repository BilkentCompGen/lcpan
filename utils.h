#ifndef __UTILS_H__
#define __UTILS_H__

#include "struct_def.h"
#include "lps.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define MARGIN 4

/**
 * @brief Frees memory allocated for the opt_arg structure.
 *
 * @param args A pointer to the `opt_arg` structure to be freed.
 */
void free_opt_arg(struct opt_arg *args);

/**
 * @brief Frees memory allocated for the ref_seq structure.
 *
 * @param seqs A pointer to the `ref_seq` structure to be freed.
 */
void free_ref_seq(struct ref_seq *seqs);

/**
 * @brief Prints reference sequences and their LCP cores in rGFA format.
 *
 * This function outputs the reference sequences and their Locally Consistent
 * Parsing (LCP) cores in rGFA/GFA format. The cores are printed as segments and
 * links, depending on the provided options.
 *
 * @param seqs     A pointer to the `ref_seq` structure containing the reference
 *                 sequences and their processed LCP cores.
 * @param is_rgfa  A flag indicating whether to print in rGFA format (1 for rGFA,
 *                 0 for GFA).
 * @param out      A file pointer to the output file where the formatted segments
 *                 and links will be written.
 */
void print_ref_seq(struct ref_seq *seqs, int is_rgfa, FILE *out);

/**
 * @brief Simulates an alternate sequence, extracts and refines LCP cores, and prints bubbles.
 *
 * This function simulates the alternate sequence based on the original sequence
 * and alternate tokens. It extracts LCP cores, refines them with the original
 * sequence, and prints the remaining differences as bubbles in rGFA format.
 * Results and errors are logged to the specified output files.
 *
 * @param chrom            A pointer to the `chr` structure representing the chromosome.
 * @param org_seq          The original sequence string.
 * @param alt_token        The alternate sequence token.
 * @param start_loc        The starting location of the variation.
 * @param lcp_level        The LCP level for refining cores.
 * @param core_id_index    A pointer to track the current core ID index.
 * @param failed_var_count A pointer to an integer to count failed variations.
 * @param bubble_count     A pointer to an integer to count bubble structures.
 * @param is_rgfa          A flag indicating if the output should be in rGFA format.
 * @param out              A file pointer for the output file to write results.
 * @param out_err          A file pointer for the error file to log issues.
 */
void variate(struct chr *chrom, const char *org_seq, const char *alt_token, uint64_t start_loc, int lcp_level, uint64_t *core_id_index, int* failed_var_count, int* bubble_count, int is_rgfa, FILE *out, FILE *out_err);

#endif