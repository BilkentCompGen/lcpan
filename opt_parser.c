#include "opt_parser.h"

void validate_file(const char *filename) {
	FILE *file = fopen(filename, "r");
	if (file == NULL) {
		fprintf(stderr, "Couldn't open %s\n", filename); 
		exit(EXIT_FAILURE);
	}
	fclose(file);
}

void parse_opts(int argc, char* argv[], struct opt_arg *args) {
	int opt;
    args->core_id_index = 0;
    args->lcp_level = DEFAULT_LCP_LEVEL;
    args->bubble_count = 0;

    int long_index;
    struct option long_options[] = {
        {"fasta", required_argument, NULL, 1},
        {"vcf", required_argument, NULL, 2},
        {"output", required_argument, NULL, 3},
        {"level", required_argument, NULL, 4},
        {NULL, 0, NULL, 0}
    };

    while ((opt = getopt_long(argc, argv, "f:v:l:o:", long_options, &long_index)) != -1) {
        switch (opt) {
		case 1:
			args->fasta_path = optarg; // reference file
            break;
        case 'f':
            args->fasta_path = optarg;
            break;
		case 2:
			args->vcf_path = optarg; // vcf file
            break;
        case 'v':
            args->vcf_path = optarg;
            break;
		case 3:
			args->rgfa_path = optarg; // output file
            break;
        case 'o':
            args->rgfa_path = optarg;
            break;
		case 4:
			args->lcp_level = atoi(optarg); // lcp level
            break;
		case 'l':
            args->lcp_level = atoi(optarg); 
            break;
        default:
            fprintf(stderr, "Error parsing options\n");
            return;
        }
    }

	validate_file(args->fasta_path);
	validate_file(args->vcf_path);

    char *fai_path = malloc(strlen(args->fasta_path) + 5);
    if (fai_path == NULL) {
        fprintf(stderr, "Memory allocation failed\n");
        return;
    }
    sprintf(fai_path, "%s.fai", args->fasta_path);
    args->fasta_fai_path = fai_path;
    
    validate_file(args->fasta_fai_path);
}