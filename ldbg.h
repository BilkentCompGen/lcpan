#ifndef __LDBG_H__
#define __LDBG_H__

#include "struct_def.h"
#include "utils.h"
#include "tpool.h"
#include <stdio.h>
#include <string.h>

/**
 * @brief Prints reference sequences and their LCP cores in rGFA format.
 *
 * This function outputs the reference sequences (LDBG) and their Locally Consistent
 * Parsing (LCP) cores in rGFA/GFA format. The cores are printed as segments and
 * links, depending on the provided options.
 *
 * @param seqs       A pointer to the `ref_seq` structure containing the reference
 *                   sequences and their processed LCP cores.
 * @param out        A file pointer to the output file where the formatted segments
 *                   and links will be written.
 */
void ldbg_print_ref_seq(struct ref_seq *seqs, FILE *out);

#endif