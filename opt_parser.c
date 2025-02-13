#include "opt_parser.h"

void validate_file(const char *filename, const char* type) {
	FILE *file = fopen(filename, "r");
	if (file == NULL) {
		fprintf(stderr, "[ERROR] Couldn't open %s (%s)\n", type, filename); 
		exit(EXIT_FAILURE);
	}
	fclose(file);
}

int summarize(struct opt_arg *args) {
    printf("[INFO] Ref: %s\n", args->fasta_path);
    if (args->program == VG) {
        printf("[INFO] VCF: %s\n", args->vcf_path);
    }
    printf("[INFO] Output: %s\n", args->gfa_path);
    printf("[INFO] GFA: %s\n", args->is_rgfa ? "rGFA" : "GFA");
    printf("[INFO] Non-Overlapping: %s\n", args->no_overlap ? "yes" : "no");
    printf("[INFO] LCP level: %d\n", args->lcp_level);
    printf("[INFO] Thread number: %d\n", args->thread_number);
    if (args->prefix != NULL) {
        printf("[INFO] Prefix: %s\n", args->prefix);
    } else {
        printf("[INFO] Prefix: lcpan\n");
    }
    return 1;
}

void printUsage() {
    fprintf(stderr, "Usage: ./lcpan [PROGRAM] [OPTIONS]\n\n");
    fprintf(stderr, "[PROGRAM]: \n");
    fprintf(stderr, "\t-vg:         Uses a variation graph-based approach.\n");
    fprintf(stderr, "\t-ldbg:       Uses LCP based de-Bruijn graph approach in construction.\n");
    // fprintf(stderr, "\t-aloe-vera:  Uses progressive genome alignment.\n");
}

void parse_opts(int argc, char* argv[], struct opt_arg *args) {

    if (argc<2) {
        printUsage();
        exit(EXIT_FAILURE);
    }

    if (strcmp(argv[1], "-vg") == 0) {
        if (argc<7) {
            fprintf(stderr, "Format: ./lcpan -vg -f ref.fa -v var.vcf -o out.rgfa [OPTIONS]\n");
            exit(EXIT_FAILURE);
        }
        args->program = VG;
    } else if (strcmp(argv[1], "-ldbg") == 0) {
        if (argc<5) {
            fprintf(stderr, "Format: ./lcpan -ldbg -f ref.fa -o out.rgfa [OPTIONS]\n");
            exit(EXIT_FAILURE);
        }
        args->program = LDBG;
    } 
    // else if (strcmp(argv[1], "-aloe-vera") == 0) {
    //     if (argc<5) {
    //         fprintf(stderr, "Format: ./lcpan -aloe-vera -f ref.fa -o out.rgfa [OPTIONS]\n");
    //         exit(EXIT_FAILURE);
    //     }
    //     args->program = ALOEVERA;
    // } 
    else {
        printUsage();
        exit(EXIT_FAILURE);
    }

    optind = 2;
    
	int opt;
    args->core_id_index = 0;
    args->lcp_level = DEFAULT_LCP_LEVEL;
    args->thread_number = DEFAULT_THREAD_NUMBER;
    args->failed_var_count = 0;
    args->invalid_line_count = 0;
    args->bubble_count = 0;
    args->is_rgfa = 1;
    args->no_overlap = 1;
    args->prefix = NULL;
    args->tload_factor = THREAD_POOL_FACTOR;
    args->verbose = 0;

    int long_index;
    struct option long_options[] = {
        {"ref", required_argument, NULL, 1},
        {"vcf", required_argument, NULL, 2},
        {"output", required_argument, NULL, 3},
        {"level", required_argument, NULL, 4},
        {"thread", required_argument, NULL, 5},
        {"rgfa", no_argument, NULL, 6},
        {"gfa", no_argument, NULL, 7},
        {"prefix", required_argument, NULL, 8},
        {"tload-factor", required_argument, NULL, 9},
        {"verbose", no_argument, NULL, 10},
        {NULL, 0, NULL, 0}
    };

    while ((opt = getopt_long(argc, argv, "r:v:o:l:t:p:s", long_options, &long_index)) != -1) {
        switch (opt) {
		case 1:
			args->fasta_path = optarg; // reference file
            break;
        case 'r':
            args->fasta_path = optarg;
            break;
		case 2:
			args->vcf_path = optarg; // vcf file
            break;
        case 'v':
            args->vcf_path = optarg;
            break;
		case 3:
			args->gfa_path = optarg; // output file
            break;
        case 'o':
            args->gfa_path = optarg;
            break;
		case 4:
			args->lcp_level = atoi(optarg); // lcp level
            break;
		case 'l':
            args->lcp_level = atoi(optarg); 
            break;
        case 5:
			args->thread_number = atoi(optarg); // thread number
            break;
        case 't':
            args->thread_number = atoi(optarg);
            break;
        case 6:
            args->is_rgfa = 1;
            break;
        case 7:
            args->is_rgfa = 0;
            break;
        case 'p':
            args->prefix = optarg;
            break;
        case 8:
            args->prefix = optarg;
            break;
        case 's':
            args->no_overlap = 0;
            break;
        case 9:
            args->tload_factor = atoi(optarg);
            break;
        case 10:
            args->verbose = 1;
            break;
        default:
            exit(EXIT_FAILURE);
        }
    }

    if (args->fasta_path == NULL) {
        fprintf(stderr, "[ERROR] Missing reference file.\n");
        exit(EXIT_FAILURE);
    }
    if (args->program == VG && args->vcf_path == NULL) {
        fprintf(stderr, "[ERROR] Missing VCF file.\n");
        exit(EXIT_FAILURE);
    }

	validate_file(args->fasta_path, "fa");
    if (args->program == VG) {
        validate_file(args->vcf_path, "vcf");
    }

    char *fai_path = malloc(strlen(args->fasta_path) + 5);
    if (fai_path == NULL) {
        fprintf(stderr, "[ERROR] Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    sprintf(fai_path, "%s.fai", args->fasta_path);
    args->fasta_fai_path = fai_path;
    
    validate_file(args->fasta_fai_path, "fai");

    (void)(args->verbose && summarize(args));
}