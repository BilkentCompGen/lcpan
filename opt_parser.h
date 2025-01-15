#ifndef __OPT_PARSER_H__
#define __OPT_PARSER_H__

#include "struct_def.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>

#define DEFAULT_LCP_LEVEL 4

void parse_opts(int argc, char* argv[], struct opt_arg *args);

#endif


