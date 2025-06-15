#include "opt_parser.h"

int ends_with(const char *str, const char *suffix) {
    size_t str_len = strlen(str);
    size_t suffix_len = strlen(suffix);

    if (str_len < suffix_len) {
        return 0;
    }

    return strcmp(str + str_len - suffix_len, suffix) == 0;
}

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
    if (args->program == VG || args->program == VGX) {
        printf("[INFO] VCF: %s\n", args->vcf_path);
    }
    printf("[INFO] Output: %s\n", args->gfa_path);
    printf("[INFO] GFA: %s, NonOv/Ov: %s, LCP level: %d, thd: %d\n", args->is_rgfa ? "rGFA" : "GFA", args->no_overlap ? "NonOv" : "Ov", args->lcp_level, args->thread_number);
    return 1;
}

void printOptions(void) {
    fprintf(stderr, "[Options]:\n");
    fprintf(stderr, "\t--ref | -r          Reference FASTA File. (.fai should be present)\n");
    fprintf(stderr, "\t--vcf | -v          VCF File.\n");
    fprintf(stderr, "\t--prefix | -p       Prefix to the log and output files. [Default: lcpan]\n");
    fprintf(stderr, "\t--level | -l        LCP Level. [Default: %d]\n", DEFAULT_LCP_LEVEL);
    fprintf(stderr, "\t--thread | -t       Thread Number. [Default: %d]\n", DEFAULT_THREAD_NUMBER);
    fprintf(stderr, "\t--rgfa | --gfa      Output Format. [Default: rGFA]\n");
    fprintf(stderr, "\t--no-overlap | -s   Allow Overlap. [Default: No]\n");
    fprintf(stderr, "\t--skip-masked       Skip Masked Chars (N). [Default: No]\n");
    fprintf(stderr, "\t--tload-factor      Number of elements that can be stored at pool at once. [Defautl: %d]\n", THREAD_POOL_FACTOR);
    fprintf(stderr, "\t--verbose  Verbose  [Default: false]\n");
}

void printUsage(void) {
    fprintf(stderr, "Usage: ./lcpan [PROGRAM] [OPTIONS]\n\n");
    fprintf(stderr, "[PROGRAM]: \n");
    fprintf(stderr, "\t-vg:         Uses a variation graph-based approach.\n");
    fprintf(stderr, "\t-vgx:        Uses a expanded variation graph-based approach.\n");
    fprintf(stderr, "\t-lbdg:       Uses LCP based de-Bruijn graph approach in construction.\n");
    // fprintf(stderr, "\t-aloe-vera:  Uses progressive genome alignment.\n");
}

void free_opt_arg(struct opt_arg *args) {
    free(args->fasta_fai_path);
    free(args->gfa_path);
	// the rest of the args char * will be freed by getops. hence, no need to free them
}

void parse_opts(int argc, char* argv[], struct opt_arg *args) {

    if (argc<2) {
        printUsage();
        exit(EXIT_FAILURE);
    }

    if (strcmp(argv[1], "-vg") == 0) {
        if (argc<6) {
            fprintf(stderr, "Format: ./lcpan -vg -r ref.fa -v var.vcf [OPTIONS]\n");
            exit(EXIT_FAILURE);
        }
        args->program = VG;
    }
    else if (strcmp(argv[1], "-vgx") == 0) {
        if (argc<6) {
            fprintf(stderr, "Format: ./lcpan -vgx -r ref.fa -v var.vcf [OPTIONS]\n");
            exit(EXIT_FAILURE);
        }
        args->program = VGX;
    } else if (strcmp(argv[1], "-lbdg") == 0) {
        if (argc<4) {
            fprintf(stderr, "Format: ./lcpan -lbdg -r ref.fa [OPTIONS]\n");
            exit(EXIT_FAILURE);
        }
        args->program = LBDG;
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
    args->core_id_index = 1;
    args->lcp_level = DEFAULT_LCP_LEVEL;
    args->thread_number = DEFAULT_THREAD_NUMBER;
    args->failed_var_count = 0;
    args->invalid_line_count = 0;
    args->bubble_count = 0;
    args->is_rgfa = 1;
    args->no_overlap = 1;
    args->skip_masked = 0;
    args->prefix = NULL;
    args->tload_factor = THREAD_POOL_FACTOR;
    args->verbose = 0;

    int long_index;
    struct option long_options[] = {
        {"ref", required_argument, NULL, 1},
        {"vcf", required_argument, NULL, 2},
        {"prefix", required_argument, NULL, 3},
        {"level", required_argument, NULL, 4},
        {"thread", required_argument, NULL, 5},
        {"rgfa", no_argument, NULL, 6},
        {"gfa", no_argument, NULL, 7},
        {"skip-masked", no_argument, NULL, 8},
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
			args->prefix = optarg; // prefix
            break;
        case 'p':
            args->prefix = optarg;
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
        case 's':
            args->no_overlap = 0;
            break;
        case 8:
            args->skip_masked = 1;
            break;
        case 9:
            args->tload_factor = atoi(optarg);
            break;
        case 10:
            args->verbose = 1;
            break;
        default:
            fprintf(stderr, "[ERROR] Invalid option %c\n", opt);
            printOptions();
            exit(EXIT_FAILURE);
        }
    }

    if (!args->is_rgfa) {
        args->skip_masked = 0;
    }

    if (args->fasta_path == NULL) {
        fprintf(stderr, "[ERROR] Missing reference file.\n");
        exit(EXIT_FAILURE);
    }
    if ((args->program == VG || args->program == VGX) && args->vcf_path == NULL) {
        fprintf(stderr, "[ERROR] Missing VCF file.\n");
        exit(EXIT_FAILURE);
    }

	validate_file(args->fasta_path, "fa");
    if (args->program == VG || args->program == VGX) {
        if (args->program == VG) {
            args->no_overlap = 1;
        }
        validate_file(args->vcf_path, "vcf");
    }

    char *fai_path = (char *) malloc(strlen(args->fasta_path)+5);
    if (fai_path == NULL) {
        fprintf(stderr, "[ERROR] Memory allocation failed\n");
        exit(EXIT_FAILURE);
    }
    sprintf(fai_path, "%s.fai", args->fasta_path);
    args->fasta_fai_path = fai_path;
    
    validate_file(args->fasta_fai_path, "fai");

    if (args->prefix == NULL) {
        args->gfa_path = strdup(args->is_rgfa ? "lcpan.rgfa" : "lcpan.gfa");
        if (!args->gfa_path) {
            fprintf(stderr, "[ERROR] strdup failed");
            exit(EXIT_FAILURE);
        }
    } else {
        if (args->is_rgfa) {
            args->gfa_path = (char *) malloc(strlen(args->prefix)+6);
            if (!args->gfa_path) {
                fprintf(stderr, "[ERROR] malloc failed");
                exit(EXIT_FAILURE);
            }
            snprintf(args->gfa_path, strlen(args->prefix)+6, "%s.rgfa", args->prefix);
        } else {
            args->gfa_path = (char *) malloc(strlen(args->prefix)+5);
            if (!args->gfa_path) {
                fprintf(stderr, "[ERROR] malloc failed");
                exit(EXIT_FAILURE);
            }
            snprintf(args->gfa_path, strlen(args->prefix)+5, "%s.gfa", args->prefix);
        }
    }

    if (args->is_rgfa != ends_with(args->gfa_path, ".rgfa")) {
        fprintf(stderr, "[WARN] Output format is %s but output file is %s\n", args->is_rgfa ? "rGFA" : "GFA", args->gfa_path);
    }

    (void)(args->verbose && summarize(args));
}
