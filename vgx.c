#include "vgx.h"

void line_queue_init(struct line_queue *queue, int capacity) {
    queue->lines = (char **)malloc(capacity * sizeof(char *));
    queue->size = 0;
    queue->capacity = capacity;
    queue->front = 0;
    queue->rear = 0;
}

void line_queue_push(struct line_queue *queue, char *line, pthread_mutex_t *queue_mutex, pthread_cond_t *cond_not_full, pthread_cond_t *cond_not_empty) {
    pthread_mutex_lock(queue_mutex);
    while (queue->size == queue->capacity) {
        pthread_cond_wait(cond_not_full, queue_mutex);
    }
    queue->lines[queue->rear] = line;
    queue->rear = (queue->rear + 1) % queue->capacity;
    queue->size++;
    pthread_cond_signal(cond_not_empty);
    pthread_mutex_unlock(queue_mutex);
}

char *line_queue_pop(struct line_queue *queue, pthread_mutex_t *queue_mutex, pthread_cond_t *cond_not_full, pthread_cond_t *cond_not_empty, int *exit_signal) {
    pthread_mutex_lock(queue_mutex);
    while (queue->size == 0 && *(exit_signal) == 0) {
        pthread_cond_wait(cond_not_empty, queue_mutex);
    }
    if (*(exit_signal) == 1) {
        pthread_mutex_unlock(queue_mutex);
        return NULL;
    }
    char *line = queue->lines[queue->front];
    queue->front = (queue->front + 1) % queue->capacity;
    queue->size--;
    pthread_cond_signal(cond_not_full);
    pthread_mutex_unlock(queue_mutex);
    return line;
}

void vgx_variate_snp(struct t_arg *t_args, const struct chr *chrom, const char *alt_token, const char *seq_name, int order, uint64_t start_loc, uint64_t end_loc, uint64_t splitting_core_id, uint64_t merging_core_id, uint64_t marginal_start, uint64_t marginal_end, int merge_overlap) {

	// print new node connecting directly from latest core before and first core after the alternating token
    print_seq3_vg(t_args->core_id_index, chrom->seq+marginal_start, start_loc-marginal_start, 
                                         alt_token, strlen(alt_token), 
                                         chrom->seq+end_loc, marginal_end-end_loc,
                                         seq_name, order, marginal_start, 1, t_args->is_rgfa, t_args->out1);
    // print splitting link
    print_link(splitting_core_id, '+', t_args->core_id_index, '+', 0, t_args->out1);
	// print merging link
    print_link(t_args->core_id_index, '+', merging_core_id, '+', merge_overlap, t_args->out1);
    t_args->core_id_index++;
}

void vgx_variate_sv(struct t_arg *t_args, const struct chr *chrom, const char *alt_token, const char *seq_name, int order, uint64_t start_loc, uint64_t end_loc, uint64_t splitting_core_id, uint64_t merging_core_id, uint64_t marginal_start, uint64_t marginal_end, int merge_overlap) {

    uint64_t alt_len = strlen(alt_token);
 
	struct lps substr;
	init_lps_offset(&substr, alt_token, alt_len, 0);
	lps_deepen(&substr, t_args->lcp_level);

    if (substr.size == 0) {
        // print first splitting node and link from reference graph	
        print_seq3_vg(t_args->core_id_index, chrom->seq+marginal_start, start_loc-marginal_start,
                                             alt_token, strlen(alt_token),
                                             chrom->seq+end_loc, marginal_end-end_loc,
                                             seq_name, order, marginal_start, 1, t_args->is_rgfa, t_args->out1);
        // print splitting link
        print_link(splitting_core_id, '+', t_args->core_id_index, '+', 0, t_args->out1);
        // print merging link
        print_link(t_args->core_id_index, '+', merging_core_id, '+', merge_overlap, t_args->out1);
        t_args->core_id_index++;
    } else {
        // print new node in between latest core before alternating core and first core in alternating core
        uint64_t prev_core_id = splitting_core_id;
        if (start_loc-marginal_start+substr.cores[0].start) {
            print_seq2_vg(t_args->core_id_index, chrom->seq+marginal_start, start_loc-marginal_start,
                                                 alt_token, substr.cores[0].start,
                                                 seq_name, order, marginal_start, 1, t_args->is_rgfa, t_args->out1);
            // print link between splitting segment with reference
            print_link(prev_core_id, '+', t_args->core_id_index, '+', 0, t_args->out1);
            prev_core_id = t_args->core_id_index;
            t_args->core_id_index++;
        }

        // print LCP cores
        if (t_args->no_overlap) {
            uint64_t prev_end = substr.cores[0].start;
            for (int i=0; i<substr.size; i++) {
                int start = maximum(prev_end, substr.cores[i].start); // if alt seq have gaps (NNN)
                int core_len = substr.cores[i].end-start;
                print_seq_vg(t_args->core_id_index, alt_token+start, core_len, seq_name, order, start_loc+start, 1, t_args->is_rgfa, t_args->out1);
                print_link(prev_core_id, '+', t_args->core_id_index, '+', 0, t_args->out1);
                prev_core_id = t_args->core_id_index;
                t_args->core_id_index++;
                prev_end = substr.cores[i].end;
            }
        } else {
            uint64_t prev_end = substr.cores[0].start;
            for (int i=0; i<substr.size; i++) {
                int start = substr.cores[i].start;
                int core_len = substr.cores[i].end-start;
                int overlap = prev_end >= substr.cores[i].start ? prev_end-substr.cores[i].start : 0;
                print_seq_vg(t_args->core_id_index, alt_token+start, core_len, seq_name, order, start_loc+start, 1, t_args->is_rgfa, t_args->out1);
                print_link(prev_core_id, '+', t_args->core_id_index, '+', overlap, t_args->out1);
                prev_core_id = t_args->core_id_index;
                t_args->core_id_index++;
                prev_end = substr.cores[i].end;
            }
        }

        if (alt_len-substr.cores[substr.size-1].end+marginal_end-end_loc) {
            // create merging segment in between last core in alternating token and reference sequence.
            print_seq2_vg(t_args->core_id_index, alt_token+substr.cores[substr.size-1].end, alt_len-substr.cores[substr.size-1].end,
                                                 chrom->seq+end_loc, marginal_end-end_loc,
                                                 seq_name, order, marginal_start, 1, t_args->is_rgfa, t_args->out1);

            print_link(prev_core_id, '+', t_args->core_id_index, '+', 0, t_args->out1);
            prev_core_id = t_args->core_id_index;
            t_args->core_id_index++;
        }

        // print merging link
        print_link(prev_core_id, '+', merging_core_id, '+', merge_overlap, t_args->out1);
    }

	free_lps(&substr);
}

uint64_t vgx_variate(struct t_arg *t_args, const struct chr *chrom, const char *org_seq, const char *alt_token, const char *seq_name, int order, uint64_t start_loc, uint64_t start_index) {

	// decide on boundaries
	uint64_t end_loc = start_loc+strlen(org_seq);
    uint64_t latest_core_index;
	uint64_t first_core_after;
	find_boundaries(start_loc, end_loc, chrom, start_index, &latest_core_index, &first_core_after);

	if (latest_core_index == 0 || chrom->cores[latest_core_index].end < chrom->cores[latest_core_index+1].start)  {
        t_args->failed_var_count += 1;
        pthread_mutex_lock(t_args->out_log_mutex);
        fprintf(t_args->out2, "VARIATE-MARGIN-START:\tCHROM: %s,\tPOSITION: %ld,\tORG: %s,\tALT: %s,\tlatest_core_index: %ld\n", chrom->seq_name, start_loc, org_seq, alt_token, latest_core_index);
        fflush(t_args->out2);
        pthread_mutex_unlock(t_args->out_log_mutex);
		return start_index;
	}
	if (first_core_after+1 >= (uint64_t)chrom->cores_size || chrom->cores[first_core_after-1].end < chrom->cores[first_core_after].start) {
        t_args->failed_var_count += 1;
        pthread_mutex_lock(t_args->out_log_mutex);
        fprintf(t_args->out2, "VARIATE-MARGIN-END:\tCHROM: %s,\tPOSITION: %ld,\tORG: %s,\tALT: %s,\tfirst_core_after: %ld,\tcores_size: %d\n", chrom->seq_name, start_loc, org_seq, alt_token, first_core_after, chrom->cores_size);
        fflush(t_args->out2);
        pthread_mutex_unlock(t_args->out_log_mutex);
		return start_index;
	}

    // get the chromosome
    uint64_t marginal_start = chrom->cores[latest_core_index].end; // start of variation
    uint64_t marginal_end = chrom->cores[first_core_after-1].end;  // end of variation
    uint64_t splitting_core_id = chrom->cores[latest_core_index].id;
    uint64_t merging_core_id = chrom->cores[first_core_after].id;
    int merge_overlap = chrom->cores[first_core_after-1].end - chrom->cores[first_core_after].start;

    t_args->bubble_count += 1;

    if (end_loc-start_loc<SV_LEN_BOUNDARY && strlen(alt_token)<SV_LEN_BOUNDARY) {
        vgx_variate_snp(t_args, chrom, alt_token, seq_name, order, start_loc, end_loc, splitting_core_id, merging_core_id, marginal_start, marginal_end, merge_overlap);
    } else {
        vgx_variate_sv(t_args, chrom, alt_token, seq_name, order, start_loc, end_loc, splitting_core_id, merging_core_id, marginal_start, marginal_end, merge_overlap);
    }

    return latest_core_index;
}

void vgx_read_vcf_thd(void *args) {

    struct t_arg *t_args = (struct t_arg *)args;
    struct line_queue *queue = (struct line_queue *)t_args->queue;

    uint64_t latest_core_index = 0;
    int latest_chrom_index = 0;

    while (1) {
        // split the line by tab characters
        char *line = line_queue_pop(queue, t_args->queue_mutex, t_args->cond_not_full, t_args->cond_not_empty, t_args->exit_signal);
        if (line == NULL) {
            break;
        }

        char *chrom, *index, *id, *seq, *alt;
        long offset;

        char *saveptr;
        chrom = strtok_r(line, "\t", &saveptr); // get chromosome name
        index = strtok_r(NULL, "\t", &saveptr); // get index  
        if (index == NULL) {
            t_args->invalid_line_count += 1;
            pthread_mutex_lock(t_args->out_log_mutex);
            fprintf(t_args->out2, "VCF: no index at line: %s\n", line);
            fflush(t_args->out2);
            pthread_mutex_unlock(t_args->out_log_mutex);
            free(line);
            continue;
        }

        offset = strtol(index, NULL, 10) - 1; // get offset
        id = strtok_r(NULL, "\t", &saveptr);  // skip ID
        seq = strtok_r(NULL, "\t", &saveptr); // get sequence
        alt = strtok_r(NULL, "\t", &saveptr); // get ALT alleles

        int chrom_index = -1;
        for (int i=0; i<t_args->seqs->size; i++) {
            if (strcmp(chrom, t_args->seqs->chrs[i].seq_name) == 0) {
                chrom_index = i;
                break;
            }
        }

        if (chrom_index == -1) {
            t_args->invalid_line_count += 1;
            pthread_mutex_lock(t_args->out_log_mutex);
            fprintf(t_args->out2, "VCF: Couldn't locate chrom %s from VCF in reference\n", chrom);
            fflush(t_args->out2);
            pthread_mutex_unlock(t_args->out_log_mutex);
            free(line);
            continue;
        }

        if (latest_chrom_index != chrom_index) {
            latest_chrom_index = chrom_index;
            latest_core_index = 0;
        }

        char *alt_saveptr;
        char *alt_token = strtok_r(alt, ",", &alt_saveptr); // split ALT alleles by comma
        int order = 0;

        while (alt_token != NULL) {
            char *alt_token_copy = strdup(alt_token); 
            if (alt_token_copy == NULL) {
                t_args->failed_var_count += 1;
                pthread_mutex_lock(t_args->out_log_mutex);
                fprintf(t_args->out2, "VCF: Memory allocation failed for alt_token_copy.\n");
                fflush(t_args->out2);
                pthread_mutex_unlock(t_args->out_log_mutex);
                continue;
            }
            latest_core_index = vgx_variate(t_args, &(t_args->seqs->chrs[chrom_index]), seq, alt_token_copy, id, order, offset, latest_core_index);
            free(alt_token_copy);

            alt_token = strtok_r(NULL, ",", &alt_saveptr);
            order++;
        }

        free(line);
    }
}

void vgx_read_vcf(struct opt_arg *args, struct ref_seq *seqs) {

    printf("[INFO] Processing variations...\n");

    FILE *out_log;
    if (args->prefix == NULL) {
        char out_err_filename[10];
        snprintf(out_err_filename, sizeof(out_err_filename), "lcpan.log");
        out_log = fopen(out_err_filename, "w");
    } else {
        char out_err_filename[strlen(args->prefix)+5];
        snprintf(out_err_filename, sizeof(out_err_filename), "%s.log", args->prefix);
        out_log = fopen(out_err_filename, "w");
    }
    
    if (out_log == NULL) {
        fprintf(stderr, "Couldn't open error log file\n");
        exit(EXIT_FAILURE);
    }
    fprintf(out_log, "vgx\n");
    fprintf(out_log, "%s\n", args->gfa_path);
    fprintf(out_log, "%d\n", args->thread_number);

    struct t_arg *t_args = (struct t_arg*)malloc(args->thread_number * sizeof(struct t_arg));
    pthread_mutex_t queue_mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t out_log_mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_cond_t cond_not_full = PTHREAD_COND_INITIALIZER;
    pthread_cond_t cond_not_empty = PTHREAD_COND_INITIALIZER;
    int exit_signal = 0;

    pthread_mutex_init(&queue_mutex, NULL);
    pthread_mutex_init(&out_log_mutex, NULL);
    pthread_cond_init(&cond_not_full, NULL);
    pthread_cond_init(&cond_not_empty, NULL);

    struct line_queue queue;
    line_queue_init(&queue, args->tload_factor * args->thread_number);

    for (int i=0; i<args->thread_number; i++) {
        char indexed_filename[strlen(args->gfa_path)+5];
        snprintf(indexed_filename, sizeof(indexed_filename), "%s.%d", args->gfa_path, i+1);
        FILE *out_thd = fopen(indexed_filename, "w");
        if (out_thd == NULL) {
            fprintf(stderr, "Couldn't open output file %s\n", indexed_filename);
            exit(EXIT_FAILURE);
        }

        t_args[i].core_id_index = ((uint64_t)(i)+1) << 32;
        t_args[i].thread_id = i+1;
        t_args[i].lcp_level = args->lcp_level;
        t_args[i].is_rgfa = args->is_rgfa;
        t_args[i].no_overlap = args->no_overlap;
        t_args[i].seqs = seqs;
        t_args[i].failed_var_count = 0;
        t_args[i].invalid_line_count = 0;
        t_args[i].bubble_count = 0;
        t_args[i].out1 = out_thd;
        t_args[i].out2 = out_log;
        t_args[i].queue = (void *)&(queue);
        t_args[i].queue_mutex = &queue_mutex;
        t_args[i].out_log_mutex = &out_log_mutex;
        t_args[i].cond_not_full = &cond_not_full;
        t_args[i].cond_not_empty = &cond_not_empty;
        t_args[i].exit_signal = &exit_signal;
    }

    struct tpool *tm;
    tm = tpool_create(args->thread_number);

    for (int i=0; i<args->thread_number; i++) {
        tpool_add_work(tm, vgx_read_vcf_thd, t_args+i);
    }

    uint64_t current_size = 1048576;
    char *line = (char *)malloc(current_size);

    FILE *file = fopen(args->vcf_path, "r");
    if (file == NULL) {
        fprintf(out_log, "VCF: Couldn't open file %s\n", args->vcf_path);
        exit(EXIT_FAILURE);
    }

    int line_count = 0;
    while (fgets(line, current_size, file) != NULL) {
        size_t len = strlen(line);
        int skip_line = 0;
        
        while (len == current_size - 1 && line[len - 1] != '\n') {
            current_size *= 2;
            char *temp_line = (char *)realloc(line, current_size);
            if (!temp_line) {
                fprintf(out_log, "VCF: Memory reallocation failed.\n");
                free(line);
                fclose(file);
                return;
            }
            line = temp_line;

            if (fgets(line + len, current_size - len, file) == NULL) {
                skip_line = 1;
                break;
            }
            len = strlen(line);
        }
        
        if (skip_line || len < 2 || line[0] == '#')
            continue;

        if (line[len - 1] == '\n') {
            line[len - 1] = '\0';
            len--;
        }

        line_count++;
        char *queue_line = strdup(line);
        if (!queue_line) {
            args->invalid_line_count += 1;
            fprintf(out_log, "VCF: Memory allocation failed.\n");
            continue;
        }
        line_queue_push(&queue, queue_line, &queue_mutex, &cond_not_full, &cond_not_empty);
    }

    fclose(file);

    pthread_mutex_lock(&queue_mutex);
    while (queue.size > 0) { 
        pthread_cond_wait(&cond_not_full, &queue_mutex);
    }
    pthread_mutex_unlock(&queue_mutex);

    exit_signal = 1;
    pthread_cond_broadcast(&cond_not_empty);
    tpool_wait(tm);
    tpool_destroy(tm);
    pthread_mutex_destroy(&queue_mutex);
    pthread_mutex_destroy(&out_log_mutex);
    pthread_cond_destroy(&cond_not_full);
    pthread_cond_destroy(&cond_not_empty);

    for (int i=0; i<args->thread_number; i++) {
        args->failed_var_count += t_args[i].failed_var_count;
        args->invalid_line_count += t_args[i].invalid_line_count;
        args->bubble_count += t_args[i].bubble_count;
        fclose(t_args[i].out1);
    }
    free(t_args);

    while (queue.size) {
        char *line = queue.lines[queue.front];
        queue.front = (queue.front + 1) % queue.capacity;
        queue.size--;
        free(line);
    }
    free(queue.lines);
    free(line);

    fclose(out_log);

    printf("[INFO] Ended processing %d lines. \n", line_count);
}