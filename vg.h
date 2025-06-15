#ifndef __VG_H__
#define __VG_H__

#include "struct_def.h"
#include "utils.h"
#include "tpool.h"
#include <stdio.h>
#include <string.h>

#define DEFAULT_ARRAY_CAPACITY 10

/**
 * @brief Reads a VCF file, processes variations, and logs output to files.
 *
 * This function reads a VCF file specified in the provided options, processes
 * the variations for the given reference sequences, and logs the output in
 * rGFA/GFA format to the specified output file. Errors and invalid lines are logged
 * to a separate error file. Additionally, it tracks the number of failed
 * variations, invalid lines, and bubble structures encountered.
 *
 * @param args A pointer to the `opt_arg` structure containing input options, 
 *             including the path to the VCF file.
 * @param seqs A pointer to the `ref_seq` structure holding the reference 
 *             sequences for variation processing.
 */
void vg_read_vcf(struct opt_arg *args, struct ref_seq *seqs);

#endif

