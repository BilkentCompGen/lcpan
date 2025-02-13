#include "vcf_parser.h"

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

void t_read_vcf(void *args) {

    struct t_arg *t_args = (struct t_arg *)args;
    struct line_queue *queue = t_args->queue;

    while (1) {
        // split the line by tab characters
        char *line = line_queue_pop(queue, t_args->queue_mutex, t_args->cond_not_full, t_args->cond_not_empty, t_args->exit_signal);
        if (line == NULL) {
            break;
        }

        char *chrom, *index, *seq, *alt;
        long offset;

        char *saveptr;
        chrom = strtok_r(line, "\t", &saveptr); // get chromosome name
        index = strtok_r(NULL, "\t", &saveptr); // get index  
        if (index == NULL) {
            t_args->invalid_line_count += 1;
            pthread_mutex_lock(t_args->out_log_mutex);
            fprintf(t_args->out_log, "VCF: no index at line: %s\n", line);
            fflush(t_args->out_log);
            pthread_mutex_unlock(t_args->out_log_mutex);
            free(line);
            continue;
        }
        
        offset = strtol(index, NULL, 10) - 1; // get offset
        strtok_r(NULL, "\t", &saveptr);       // skip ID
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
            fprintf(t_args->out_log, "VCF: Couldn't locate chrom %s from VCF in reference\n", chrom);
            fflush(t_args->out_log);
            pthread_mutex_unlock(t_args->out_log_mutex);
            free(line);
            continue;
        }

        char *alt_saveptr;
        char *alt_token = strtok_r(alt, ",", &alt_saveptr); // split ALT alleles by comma

        while (alt_token != NULL) {
            char *alt_token_copy = strdup(alt_token); 
            if (alt_token_copy == NULL) {
                t_args->failed_var_count += 1;
                pthread_mutex_lock(t_args->out_log_mutex);
                fprintf(t_args->out_log, "VCF: Memory allocation failed for alt_token_copy.\n");
                fflush(t_args->out_log);
                pthread_mutex_unlock(t_args->out_log_mutex);
                continue;
            }
            variate(t_args, &(t_args->seqs->chrs[chrom_index]), seq, alt_token_copy, offset);
            free(alt_token_copy);

            alt_token = strtok_r(NULL, ",", &alt_saveptr);
        }

        free(line);
    }
}

void read_vcf(struct opt_arg *args, struct ref_seq *seqs) {

    printf("[INFO] Started processing variation file.\n");

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
        FILE *out = fopen(indexed_filename, "w");
        if (out == NULL) {
            fprintf(stderr, "Couldn't open output file %s\n", indexed_filename);
            exit(EXIT_FAILURE);
        }

        t_args[i].core_id_index = 0;
        t_args[i].thread_id = i+1;
        t_args[i].lcp_level = args->lcp_level;
        t_args[i].is_rgfa = args->is_rgfa;
        t_args[i].no_overlap = args->no_overlap;
        t_args[i].seqs = seqs;
        t_args[i].failed_var_count = 0;
        t_args[i].invalid_line_count = 0;
        t_args[i].bubble_count = 0;
        t_args[i].out = out;
        t_args[i].out_log = out_log;
        t_args[i].queue = &(queue);
        t_args[i].queue_mutex = &queue_mutex;
        t_args[i].out_log_mutex = &out_log_mutex;
        t_args[i].cond_not_full = &cond_not_full;
        t_args[i].cond_not_empty = &cond_not_empty;
        t_args[i].exit_signal = &exit_signal;
    }

    struct tpool *tm;

    tm = tpool_create(args->thread_number);

    for (int i=0; i<args->thread_number; i++) {
        tpool_add_work(tm, t_read_vcf, t_args+i);
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
        line_count++;
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

        char *queue_line = strdup(line);
        if (!queue_line) {
            args->invalid_line_count += 1;
            fprintf(out_log, "VCF: Memory allocation failed.\n");
            continue;
        }
        line_queue_push(&queue, queue_line, &queue_mutex, &cond_not_full, &cond_not_empty);
    }

    fclose(file);

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
        fclose(t_args[i].out);
    }

    free(queue.lines);
    free(t_args);

    fclose(out_log);

    printf("[INFO] Ended processing variation file. line_count: %d\n", line_count);
}