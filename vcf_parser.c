#include "vcf_parser.h"

void line_queue_init(struct line_queue *queue, int capacity) {
    queue->lines = (char **)malloc(capacity * sizeof(char *));
    queue->size = 0;
    queue->capacity = capacity;
    queue->front = 0;
    queue->rear = 0;
}

void line_queue_push(struct line_queue *queue, char *line, pthread_mutex_t *mutex, pthread_cond_t *cond_not_full, pthread_cond_t *cond_not_empty) {
    pthread_mutex_lock(mutex);
    while (queue->size == queue->capacity) {
        pthread_cond_wait(cond_not_full, mutex);
    }
    queue->lines[queue->rear] = line;
    queue->rear = (queue->rear + 1) % queue->capacity;
    queue->size++;
    pthread_cond_signal(cond_not_empty);
    pthread_mutex_unlock(mutex);
}

char *line_queue_pop(struct line_queue *queue, pthread_mutex_t *mutex, pthread_cond_t *cond_not_full, pthread_cond_t *cond_not_empty, int *exit_signal) {
    pthread_mutex_lock(mutex);
    while (queue->size == 0 && *(exit_signal) == 0) {
        pthread_cond_wait(cond_not_empty, mutex);
    }
    if (*(exit_signal) == 1) {
        pthread_mutex_unlock(mutex);
        return NULL;
    }
    char *line = queue->lines[queue->front];
    queue->front = (queue->front + 1) % queue->capacity;
    queue->size--;
    pthread_cond_signal(cond_not_full);
    pthread_mutex_unlock(mutex);
    return line;
}

void t_read_vcf(void *args) {

    struct t_arg *t_args = (struct t_arg *)args;

    FILE *out = t_args->out;
    FILE *out_err = t_args->out_err;
    struct line_queue *queue = t_args->queue;

    while (1) {
        // split the line by tab characters
        char *line = line_queue_pop(queue, t_args->mutex, t_args->cond_not_full, t_args->cond_not_empty, t_args->exit_signal);
        if (*(t_args->exit_signal) == 1) {
            pthread_mutex_unlock(t_args->mutex);
            break;
        }
        if (line == NULL) {
            pthread_mutex_lock(t_args->out_err_mutex);
            t_args->invalid_line_count += 1;
            fprintf(out_err, "VCF: line is empty.\n");
            pthread_mutex_unlock(t_args->out_err_mutex);
            continue;
        }

        char *chrom, *index, *seq, *alt;
        long offset;

        char *saveptr;
        chrom = strtok_r(line, "\t", &saveptr); // get chromosome name
        index = strtok_r(NULL, "\t", &saveptr); // get index  
        if (index == NULL) {
            pthread_mutex_lock(t_args->out_err_mutex);
            t_args->invalid_line_count += 1;
            fprintf(out_err, "VCF: no index at line: %s\n", line);
            pthread_mutex_unlock(t_args->out_err_mutex);
            free(line);
            continue;
        }
        
        offset = strtol(index, NULL, 10);     // get offset
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
            pthread_mutex_lock(t_args->out_err_mutex);
            t_args->invalid_line_count += 1;
            fprintf(out_err, "VCF: Couldn't locate chrom %s from VCF in reference\n", chrom);
            pthread_mutex_unlock(t_args->out_err_mutex);
            free(line);
            continue;
        }

        char *alt_saveptr;
        char *alt_token = strtok_r(alt, ",", &alt_saveptr); // split ALT alleles by comma
        
        while (alt_token != NULL) {
            variate(t_args->args, &(t_args->seqs->chrs[chrom_index]), seq, alt_token, offset, &(t_args->failed_var_count), t_args->core_id_mutex, t_args->out_err_mutex, out, out_err);
            alt_token = strtok_r(NULL, ",", &alt_saveptr);
        }

        free(line);
    }
}

void read_vcf(struct opt_arg *args, struct ref_seq *seqs, int *failed_var_count, int *invalid_line_count, FILE *out_err) {

    printf("[INFO] Started processing variation file.\n");

    struct t_arg *t_args = (struct t_arg*)malloc(args->thread_number * sizeof(struct t_arg));
    pthread_mutex_t mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t out_err_mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_mutex_t core_id_mutex = PTHREAD_MUTEX_INITIALIZER;
    pthread_cond_t cond_not_full = PTHREAD_COND_INITIALIZER;
    pthread_cond_t cond_not_empty = PTHREAD_COND_INITIALIZER;
    int exit_signal = 0;

    pthread_mutex_init(&mutex, NULL);
    pthread_mutex_init(&out_err_mutex, NULL);
    pthread_mutex_init(&core_id_mutex, NULL);
    pthread_cond_init(&cond_not_full, NULL);
    pthread_cond_init(&cond_not_empty, NULL);

    struct line_queue queue;
    line_queue_init(&queue, 2*args->thread_number);

    for (int i=0; i<args->thread_number; i++) {
        t_args[i].args = args;
        t_args[i].seqs = seqs;
        t_args[i].failed_var_count = 0;
        t_args[i].invalid_line_count = 0;
        t_args[i].out_err = out_err;
        t_args[i].queue = &(queue);
        t_args[i].mutex = &mutex;
        t_args[i].out_err_mutex = &out_err_mutex;
        t_args[i].core_id_mutex = &core_id_mutex;
        t_args[i].cond_not_full = &cond_not_full;
        t_args[i].cond_not_empty = &cond_not_empty;
        t_args[i].exit_signal = &exit_signal;

        char indexed_filename[strlen(args->gfa_path)+5];
        snprintf(indexed_filename, sizeof(indexed_filename), "%s.%d", args->gfa_path, i);
        FILE *out = fopen(indexed_filename, "w");
        if (out == NULL) {
            fprintf(stderr, "Couldn't open output file %s\n", indexed_filename);
            exit(EXIT_FAILURE);
        }
        t_args[i].out = out;
    }

    struct tpool *tm;

    tm = tpool_create(args->thread_number);

    for (int i=0; i<args->thread_number; i++) {
        tpool_add_work(tm, t_read_vcf, t_args+i);
    }

    uint64_t current_size = 1048576;
    char *line = (char *)malloc(current_size);

    FILE *file = fopen(t_args->args->vcf_path, "r");
    if (file == NULL) {
        fprintf(out_err, "VCF: Couldn't open file %s\n", t_args->args->vcf_path);
        exit(EXIT_FAILURE);
    }

    int line_count = 0;
    while (fgets(line, current_size, file) != NULL) {
        while (strlen(line) == current_size - 1 && line[current_size - 2] != '\n') {
            long current_position = ftell(file);
            fseek(file, current_position - (current_size - 1), SEEK_SET); // last char is '\0'
            
            current_size *= 2;
            line = (char *)realloc(line, current_size);
            continue; // retry to read line
        }

        line[strlen(line)-1] = '\0';
        
        if (line[0] == '#')
            continue;

        char *queue_line = (char *)strdup(line);
        line_queue_push(&queue, queue_line, &mutex, &cond_not_full, &cond_not_empty);
        line_count++;
    }

    fclose(file);

    exit_signal = 1;
    pthread_cond_broadcast(&cond_not_empty);
    tpool_wait(tm);
    tpool_destroy(tm);
    pthread_mutex_destroy(&mutex);
    pthread_mutex_destroy(&out_err_mutex);
    pthread_mutex_destroy(&core_id_mutex);
    pthread_cond_destroy(&cond_not_full);
    pthread_cond_destroy(&cond_not_empty);

    for (int i=0; i<args->thread_number; i++) {
        *failed_var_count = *failed_var_count + t_args->failed_var_count;
        *invalid_line_count = *invalid_line_count + t_args->invalid_line_count;
        fclose(t_args[i].out);
    }

    free(queue.lines);
    free(t_args);

    printf("[INFO] Ended processing variation file. line_count: %d\n", line_count);
}