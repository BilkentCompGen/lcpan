#ifndef __FASTA_PARSER__
#define __FASTA_PARSER__

#include "struct_def.h"
#include "lps.h"
#include "utils.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

/**
 * @brief Frees memory allocated for the ref_seq structure.
 *
 * @param seqs A pointer to the `ref_seq` structure to be freed.
 */
void free_ref_seq(struct ref_seq *seqs);

/**
 * @brief Reads a FASTA file and processes sequences using the LCP.
 *
 * This function reads a FASTA file specified in the provided options,
 * extracts chromosome sequences, and processes them to generate cores
 * using the Locally Consistent Parsing (LCP). The processed sequences
 * and their cores are stored in the provided `ref_seq` structure.
 *
 * @param args A pointer to the `opt_arg` structure containing input options,
 *             including the path to the FASTA file.
 * @param seqs A pointer to the `ref_seq` structure where the extracted 
 *             sequences and their LCP-processed cores will be stored.
 */
void read_fasta(struct opt_arg *args, struct ref_seq *seqs);

/**
 * @brief Prints reference sequences and their LCP cores in rGFA format.
 *
 * This function outputs the reference sequences and their Locally Consistent
 * Parsing (LCP) cores in rGFA/GFA format. The cores are printed as segments and
 * links, depending on the provided options.
 *
 * @param seqs       A pointer to the `ref_seq` structure containing the reference
 *                   sequences and their processed LCP cores.
 * @param is_rgfa    A flag indicating whether to print in rGFA format (1 for rGFA,
 *                   0 for GFA).
 * @param out        A file pointer to the output file where the formatted segments
 *                   and links will be written.
 */
void print_ref_seqs(const struct ref_seq *seqs, int is_rgfa, FILE *out);

#endif