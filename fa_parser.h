#ifndef __FASTA_PARSER__
#define __FASTA_PARSER__

#include "struct_def.h"
#include "lps.h"
#include <stdio.h>
#include <string.h>

void read_fasta(struct opt_arg *args, struct ref_seq *seqs);

#endif