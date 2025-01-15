#ifndef __VCF_PARSER_H__
#define __VCF_PARSER_H__

#include "struct_def.h"
#include "utils.h"
#include <stdio.h>
#include <string.h>

void read_vcf(struct opt_arg *args, struct ref_seq *seqs, int* failed_var_count, FILE *out, FILE *out_err);

#endif