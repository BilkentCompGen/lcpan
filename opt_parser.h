#ifndef __OPT_PARSER_H__
#define __OPT_PARSER_H__

#include "struct_def.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>

#define DEFAULT_LCP_LEVEL 4
#define DEFAULT_THREAD_NUMBER 1

/**
 * @brief Parses command-line options and arguments.
 *
 * This function processes the command-line arguments passed to the program,
 * extracts options and their associated values, and populates the provided
 * `opt_arg` structure.
 *
 * @param argc The number of arguments passed to the program.
 * @param argv The array of argument strings.
 * @param args A pointer to the `opt_arg` structure where the parsed options
 *             and their values will be stored.
 */
void parse_opts(int argc, char* argv[], struct opt_arg *args);

#endif


