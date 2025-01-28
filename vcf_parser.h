#ifndef __VCF_PARSER_H__
#define __VCF_PARSER_H__

#include "struct_def.h"
#include "utils.h"
#include "tpool.h"
#include <stdio.h>
#include <string.h>

/**
 * @brief Reads a VCF file, processes variations, and logs output to files.
 *
 * This function reads a VCF file specified in the provided options, processes
 * the variations for the given reference sequences, and logs the output in
 * rGFA/GFA format to the specified output file. Errors and invalid lines are logged
 * to a separate error file. Additionally, it tracks the number of failed
 * variations, invalid lines, and bubble structures encountered.
 *
 * @param args               A pointer to the `opt_arg` structure containing
 *                           input options, including the path to the VCF file.
 * @param seqs               A pointer to the `ref_seq` structure holding the
 *                           reference sequences for variation processing.
 * @param failed_var_count   A pointer to an integer to store the count of 
 *                           failed variations.
 * @param invalid_line_count A pointer to an integer to store the count of 
 *                           invalid lines in the VCF file.
 * @param out_err            A file pointer to the error file where errors and
 *                           invalid lines will be logged.
 */
void read_vcf(struct opt_arg *args, struct ref_seq *seqs, int *failed_var_count, int *invalid_line_count, FILE *out_err);

#endif