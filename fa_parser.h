#ifndef __FASTA_PARSER__
#define __FASTA_PARSER__

#include "struct_def.h"
#include "lps.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

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

#endif