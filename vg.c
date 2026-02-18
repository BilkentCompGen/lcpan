#include "vg.h"


// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------
//      UTILS
// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------

/**
 * Print sequence. It may happen when there is a chromosomal jump (i.e., no variation
 * on entire chromosome.)
 * 
 */
void vg_print_seq(struct chr *chrom, int is_rgfa, FILE *out_segment, FILE *out_link) {
    if (chrom->cores_size) {
        chrom->ids = NULL; // To print simple path
        const char *seq_name = chrom->seq_name;
        const char *seq = chrom->seq;
        int cores_size = chrom->cores_size;
        
        {
            const struct simple_core *temp_core = &(chrom->cores[0]);
            uint64_t start = temp_core->start;
            int seq_len = (int)(temp_core->end - start);
            print_seq(temp_core->id, seq+start, seq_len, seq_name, start, 0, is_rgfa, out_segment);
        }
        
        uint64_t prev_core_id = chrom->cores[0].id;

        for (int j=1; j<cores_size; j++) {
            const struct simple_core *temp_core = &(chrom->cores[j]);
            uint64_t temp_id = temp_core->id;
            uint64_t temp_start = temp_core->start;
            int seq_len = (int)(temp_core->end - temp_start);

            print_seq(temp_core->id, seq+temp_start, seq_len, seq_name, temp_start, 0, is_rgfa, out_segment);
            print_link(prev_core_id, '+', temp_core->id, '+', 0, out_link);
            prev_core_id = temp_id;
        }
    }
}

// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------
//      THREADS
// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------

/**
 * Naming threads for better profiling analyses.
 */
static void name_thread(const char* name) {
#if defined(__APPLE__)
    pthread_setname_np(name);

#elif defined(__linux__)
    pthread_setname_np(pthread_self(), name);

#elif defined(_WIN32)
    wchar_t wname[64];
    MultiByteToWideChar(CP_UTF8, 0, name, -1, wname, 64);
    SetThreadDescription(GetCurrentThread(), wname);

#else
    (void)name; // unsupported
#endif
}

// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------
//      HELPERS
// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------

/**
 * Set ID for the segment. Process each variation graph data (bucket) and if there is a
 * outgoing deletion from the segment's end, then assign the pre-defined ID. If not,
 * assign global/thread specific incremental ID.
 */
uint64_t set_id(vg_core_bucket_t *bucket, uint64_t segment_end, uint64_t *t_args_id) {
    for (int i = 0; i < bucket->size; i++) {
        if (bucket->items[i].start == segment_end && bucket->items[i].dir == VG_DIR_OUT && bucket->items[i].var == VG_VAR_DEL) { // if DEL start
            return bucket->items[i].id;
        }
    }
    (*t_args_id)++;
    return (*t_args_id)-1;
}

/**
 * Locate split and merge segments' IDs of the given variation. If variation starts 
 * at the beginning, then assign the split_id the previous_segment. Note that variation
 * should have start position so split_id can be assigned. Also, it should have end
 * position.
 */
void locate_ids(vg_core_bucket_t *bucket, struct simple_core *segments, int segment_count, uint64_t start, uint64_t end, uint64_t *split_id, uint64_t *merge_id) {
    int i = 0;
    *split_id = bucket->prev_id;
    *merge_id = bucket->curr_id;
    if (start != 0xFFFFFFFFFFFFFFFF) {
        for(; i < segment_count; i++) {
            if (start == segments[i].start) {
                *split_id = i ? segments[i-1].id : bucket->prev_id;
                break;
            }
        }
    }

    if (end != 0xFFFFFFFFFFFFFFFF) {
        for(; i < segment_count; i++) {
            if (end == segments[i].start) {
                *merge_id = segments[i].id;
                break;
            }
        }
    }
}

/**
 * This function adds the variations that are note part of the current vdg (i.e., outgoing).
 * Note that data is stored in sorted order, based on the positions. The element that is
 * inserted is stored as 64-bit integer, (32-bit ID + 32-bit Pos), which the ID is assigned 
 * by main thread.
 */
static inline void add_pending_var_end(uint64_t ** restrict pending_var_ends, int * restrict pending_var_ends_size, int * restrict pending_var_ends_capacity, uint64_t id, uint64_t loc) {
    int size      = *pending_var_ends_size;
    int cap       = *pending_var_ends_capacity;
    uint64_t *arr = *pending_var_ends;

    if (size == cap) {
        uint64_t *temp = (uint64_t *)realloc(arr, sizeof(uint64_t) * 2 * cap);
        if (!temp) { perror("[ERROR] Failed to increase remaining variants array."); abort(); }
        cap *= 2;
        arr = temp;
        *pending_var_ends = arr;
        *pending_var_ends_capacity = cap;
    }

    arr[size++] = (id << 32) | (uint32_t)loc;

    int i = size - 1;
    while (i > 0 && (uint32_t)arr[i - 1] > (uint32_t)arr[i]) {
        uint64_t t = arr[i - 1];
        arr[i - 1] = arr[i];
        arr[i] = t;
        i--;
    }

    *pending_var_ends_size = size;
}

/**
 * This function facilitates to move the bucket data to the next segment. Note that it 
 * does not make allocation for the bucket->items variations array. It only modifies metadata.
 * Hence, you need to make sure that the previous segment should nothing to process (split).
 */
static inline void move_to_next_core(vg_core_bucket_t *bucket, int chr_idx, int core_idx, uint64_t curr_id) {
    bucket->chr_idx   = chr_idx;
    bucket->core_idx  = core_idx;
    bucket->size      = 0;
    bucket->prev_id   = bucket->curr_id;
    bucket->curr_id   = curr_id;
}

/**
 * This function creates new bucket data. It initializes a new array to store variatons.
 * All other necessary information is also initialized.
 */
static inline vg_core_bucket_t *malloc_vg_core_bucket(int chr_idx, int core_idx, uint64_t curr_id, uint64_t prev_id) {
    vg_core_bucket_t *bucket = (vg_core_bucket_t *)malloc(sizeof(vg_core_bucket_t));
    bucket->chr_idx     = chr_idx;
    bucket->core_idx    = core_idx;
    bucket->capacity    = DEFAULT_ARRAY_CAPACITY;
    bucket->size        = 0;
    bucket->curr_id     = curr_id;
    bucket->prev_id     = prev_id;
    bucket->items       = (vg_element_t *)malloc(sizeof(vg_element_t) * bucket->capacity);
    return bucket;
}

/**
 * This function checks whether there is any space left to insert variation (element)
 * into bucket->items array. If not, it increases the array by reallocating data.
 */
static inline int check_vg_items(vg_core_bucket_t *bucket) {
    if (bucket->size == bucket->capacity) {
        bucket->capacity *= 2;
        vg_element_t *temp = (vg_element_t *)realloc(bucket->items, sizeof(vg_element_t) * bucket->capacity);
        if (temp == NULL) return 0;
        bucket->items = temp;
    }
    return 1;
}

/**
 * As for newly created data, this function add elements to bucket->items array.
 * This function basically assigns incoming variations to the segment if there is any.
 * So, threads can make necesarry linking of incoming variatons (such as del, alt...)
 */
static inline void add_element_to_bucket(vg_core_bucket_t *bucket, uint64_t *arr, int *arr_size, uint64_t start, uint64_t end) {
    int i = 0, size = *arr_size;
    // If elements in rem_arr lies in this core, add them to vg_data
    while (i < size && (uint32_t)(arr[i]) < start) {
        // this only happens when the end point of variation is in masked region (N)
        // if masked regions are represented in segments, no problem will occur
        fprintf(stderr, "[ERROR] Left element on the way. chrom: %d, var end: %u, core start: %lu\n", bucket->chr_idx, (uint32_t)(arr[i]), start);
        i++;
    }
    while (i < size && (uint32_t)(arr[i]) < end) {
        check_vg_items(bucket); // check if there is a space to add element
        bucket->items[bucket->size] = (vg_element_t){VG_DIR_INCOMING, VG_VAR_NONE, arr[i] >> 32, 0xFFFFFFFFFFFFFFFF, 0xFFFFFFFF & arr[i], NULL, NULL, 0}; // order not matter here
        bucket->size++;
        i++;
    }
    // remove first elements by shifting the data to the left
    if (i && i < size) {
        memmove(arr, arr+i, (size-i) * sizeof(uint64_t));
    }

    *arr_size = size - i;
}

static inline void vg_print_core_as_is(const struct chr *chr, int chr_idx, int core_idx, struct ref_seq *seqs, int is_rgfa, FILE *out_segment, FILE *out_link) {
    const struct simple_core *c = &chr->cores[core_idx];
    uint64_t core_id            = c->id;
    uint64_t core_start         = c->start;
    uint64_t core_end           = c->end;

    print_seq(core_id, chr->seq + core_start, core_end - core_start, chr->seq_name, core_start, 0, is_rgfa, out_segment);

    // in case it is first lcp core in chromosome
    if (core_idx) {
        print_link(chr->cores[core_idx - 1].id, '+', core_id, '+', 0, out_link);
    }

    seqs->chrs[chr_idx].ids[core_idx] = NULL;
}

// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------
//      THREAD POOL UTILS
// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------

void vg_queue_init(vg_work_queue_t *queue, int capacity) {
    queue->batch    = (vg_bucket_batch_t **)malloc(sizeof(vg_bucket_batch_t *) * capacity);
    queue->size     = 0;
    queue->capacity = capacity;
    queue->front    = 0;
    queue->rear     = 0;
}

static inline void vg_queue_push(vg_work_queue_t *queue, vg_bucket_batch_t *batch, vg_queue_sync_t *sync) {
    pthread_mutex_lock(sync->mutex);
    while (queue->size == queue->capacity) {
        pthread_cond_wait(sync->cond_not_full, sync->mutex);
    }
    queue->batch[queue->rear] = batch;
    queue->rear = (queue->rear + 1) % queue->capacity;
    queue->size++;
    pthread_cond_signal(sync->cond_not_empty);
    pthread_mutex_unlock(sync->mutex);
}

static inline vg_bucket_batch_t *vg_queue_pop(vg_work_queue_t *queue, vg_queue_sync_t *sync) {
    pthread_mutex_lock(sync->mutex);
    while (queue->size == 0 && *(sync->exit_signal) == 0) {
        pthread_cond_wait(sync->cond_not_empty, sync->mutex);
    }
    if (*(sync->exit_signal) == 1) {
        pthread_mutex_unlock(sync->mutex);
        return NULL;
    }
    vg_bucket_batch_t *bucket = queue->batch[queue->front];
    queue->front = (queue->front + 1) % queue->capacity;
    queue->size--;
    pthread_cond_signal(sync->cond_not_full);
    pthread_mutex_unlock(sync->mutex);
    return bucket;
}

// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------
//      SV HANDLER
// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------

void vg_variate_sv(struct t_arg *t_args, uint64_t split_id, uint64_t merge_id, vg_element_t *sv) {

    uint64_t alt_len = strlen(sv->seq);
 
	struct lps substr;
	init_lps_offset(&substr, sv->seq, alt_len, 0);
	lps_deepen(&substr, t_args->lcp_level);

    if (substr.size) {
        int start = sv->start;
        uint64_t prev_id = split_id;
        uint64_t prev_index = 0;
        if (substr.cores[0].start) {
            print_seq_vg(t_args->core_id_index, sv->seq, substr.cores[0].start, sv->seq_id, sv->order, start, 1, t_args->is_rgfa, t_args->out1);
            print_link(prev_id, '+', t_args->core_id_index, '+', 0, t_args->out2);
            prev_id = t_args->core_id_index;
            prev_index = substr.cores[0].start;
            t_args->core_id_index++;
        }
        int lcp_core_end_index = substr.size;
        if (substr.cores[substr.size-1].end == alt_len)
            lcp_core_end_index--;
        
        for(int i=0; i<lcp_core_end_index; i++) {
            print_seq_vg(t_args->core_id_index, sv->seq+prev_index, substr.cores[i].end-prev_index, sv->seq_id, sv->order, start+prev_index, 1, t_args->is_rgfa, t_args->out1);
            print_link(prev_id, '+', t_args->core_id_index, '+', 0, t_args->out2);
            prev_id = t_args->core_id_index;
            prev_index = substr.cores[i].end;
            t_args->core_id_index++;
        }
        print_seq_vg(sv->id, sv->seq+prev_index, alt_len-prev_index, sv->seq_id, sv->order, start+prev_index, 1, t_args->is_rgfa, t_args->out1);
        print_link(prev_id, '+', sv->id, '+', 0, t_args->out2);       
    } else {
        print_seq_vg(sv->id, sv->seq, alt_len, sv->seq_id, sv->order, sv->start, 1, t_args->is_rgfa, t_args->out1);
        print_link(split_id, '+', sv->id, '+', 0, t_args->out2);
    }

    if (merge_id) {
        print_link(sv->id, '+', merge_id, '+', 0, t_args->out2);
    }

    free(sv->seq);
    free(sv->seq_id);

	free_lps(&substr);
}

// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------
//      WORKER
// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------

void vg_read_vcf_thd(void *args) {

    struct t_arg *t_args = (struct t_arg *)args;
    vg_work_queue_t *queue = (vg_work_queue_t *)t_args->queue;

    if (t_args->out1 == NULL || t_args->out2 == NULL) {
        return;
    }

    // name thread
    char thread_name[16];
    snprintf(thread_name, sizeof(thread_name), "worker-%d", t_args->thread_id);
    name_thread(thread_name);

    time_t thread_start;
    time(&thread_start);

    while (1) {
        vg_bucket_batch_t *batch = vg_queue_pop(queue, t_args->sync);
		if (batch == NULL) break;

        for (int i = 0; i < batch->count; i++) {
            
            vg_core_bucket_t *bucket = batch->items[i];

            if (bucket->size == 0) {
                vg_print_core_as_is(&(t_args->seqs->chrs[bucket->chr_idx]), bucket->chr_idx, bucket->core_idx, t_args->seqs, t_args->is_rgfa, t_args->out1, t_args->out2);
                free(bucket->items); free(bucket);
                continue;
            }

            // split LCP core into segments
            uint64_t *split_points = (uint64_t *)malloc((2 * bucket->size + 2) * sizeof(uint64_t));
            int size = 0;
            for (int i = 0; i < bucket->size; i++) {
                if (bucket->items[i].start != 0xFFFFFFFFFFFFFFFF) split_points[size++] = bucket->items[i].start;
                if (bucket->items[i].end   != 0xFFFFFFFFFFFFFFFF) split_points[size++] = bucket->items[i].end;
            }
            split_points[size++] = t_args->seqs->chrs[bucket->chr_idx].cores[bucket->core_idx].start;
            split_points[size++] = t_args->seqs->chrs[bucket->chr_idx].cores[bucket->core_idx].end;
            quicksort(split_points, 0, size - 1); // to remove duplicate locations

            int segment_count = 0;
            for (int i = 1; i < size; i++) {
                if (split_points[segment_count] != split_points[i]) split_points[++segment_count] = split_points[i];
            }

            struct simple_core *segments = (struct simple_core *)malloc(sizeof(struct simple_core) * segment_count);
            // initialize segments data (id, start, end) and print links and segments. 
            // Note: there might be outgoing edge with pre-defined id. Make sure the ids are assigned properly
            // Note: link first segment with bucket->prev_id, give id to last element sizeof(struct simple_core)->curr_id
            if (1 < segment_count) {
                t_args->seqs->chrs[bucket->chr_idx].ids[bucket->core_idx] = (uint64_t *)malloc(sizeof(uint64_t) * segment_count);
                const char *seq                     = t_args->seqs->chrs[bucket->chr_idx].seq;
                const char *seq_name                = t_args->seqs->chrs[bucket->chr_idx].seq_name;
                const struct simple_core *curr_core = &(t_args->seqs->chrs[bucket->chr_idx].cores[bucket->core_idx]);
                uint64_t prev_segment_id            = bucket->prev_id;
                
                for (int k = 0; k < segment_count - 1; k++) {
                    uint64_t segment_id = set_id(bucket, split_points[k+1], &(t_args->core_id_index));
                    segments[k] = (struct simple_core){segment_id, split_points[k], split_points[k+1]};
                    
                    print_seq(segment_id, seq + split_points[k], split_points[k + 1] - split_points[k], seq_name, split_points[k], 0, t_args->is_rgfa, t_args->out1);
                    print_link(prev_segment_id, '+', segment_id, '+', 0, t_args->out2);
                    
                    prev_segment_id = segment_id;
                    t_args->seqs->chrs[bucket->chr_idx].ids[bucket->core_idx][k] = segment_id;
                }
                segments[segment_count - 1]       = (struct simple_core){bucket->curr_id, split_points[segment_count - 1], curr_core->end};
                
                print_seq(bucket->curr_id, seq + split_points[segment_count - 1], curr_core->end-split_points[segment_count - 1], seq_name, split_points[segment_count - 1], 0, t_args->is_rgfa, t_args->out1);
                print_link(prev_segment_id, '+', bucket->curr_id, '+', 0, t_args->out2);
                
                t_args->seqs->chrs[bucket->chr_idx].ids[bucket->core_idx][segment_count - 1] = 0;
            } else {
                const char *seq                     = t_args->seqs->chrs[bucket->chr_idx].seq;
                const char *seq_name                = t_args->seqs->chrs[bucket->chr_idx].seq_name;
                const struct simple_core *curr_core = &(t_args->seqs->chrs[bucket->chr_idx].cores[bucket->core_idx]);
                segments[0]                         = (struct simple_core){bucket->curr_id, curr_core->start, curr_core->end};
                
                print_seq(bucket->curr_id, seq+curr_core->start, curr_core->end-curr_core->start, seq_name, curr_core->start, 0, t_args->is_rgfa, t_args->out1);
                print_link(bucket->prev_id, '+', bucket->curr_id, '+', 0, t_args->out2);
                
                t_args->seqs->chrs[bucket->chr_idx].ids[bucket->core_idx] = NULL;
            }

            // iterate through every variation and print segments and links
            for (int i = 0; i < bucket->size; i++) {
                uint64_t split_id, merge_id;
                switch (bucket->items[i].dir) {
                case VG_DIR_IN:
                    if (bucket->items[i].var == VG_VAR_SNP) {
                        locate_ids(bucket, segments, segment_count, bucket->items[i].start, bucket->items[i].end, &split_id, &merge_id);
                        print_link(split_id, '+', bucket->items[i].id, '+', 0, t_args->out2);
                        print_link(bucket->items[i].id, '+', merge_id, '+', 0, t_args->out2);
                    } else if (bucket->items[i].var == VG_VAR_INS) {
                        locate_ids(bucket, segments, segment_count, bucket->items[i].start, bucket->items[i].end, &split_id, &merge_id);
                        print_link(split_id, '+', bucket->items[i].id, '+', 0, t_args->out2);
                        print_link(bucket->items[i].id, '+', merge_id, '+', 0, t_args->out2);
                    } else if (bucket->items[i].var == VG_VAR_INS_SV) {
                        locate_ids(bucket, segments, segment_count, bucket->items[i].start, bucket->items[i].end, &split_id, &merge_id);
                        vg_variate_sv(t_args, split_id, merge_id, &(bucket->items[i]));
                    } else if (bucket->items[i].var == VG_VAR_DEL) {
                        locate_ids(bucket, segments, segment_count, bucket->items[i].start, bucket->items[i].end, &split_id, &merge_id);
                        print_link(split_id, '+', merge_id, '+', 0, t_args->out2);
                    } else if (bucket->items[i].var == VG_VAR_ALT) {
                        locate_ids(bucket, segments, segment_count, bucket->items[i].start, bucket->items[i].end, &split_id, &merge_id);
                        print_link(split_id, '+', bucket->items[i].id, '+', 0, t_args->out2);
                        print_link(bucket->items[i].id, '+', merge_id, '+', 0, t_args->out2);
                    } else if (bucket->items[i].var == VG_VAR_ALT_SV) {
                        locate_ids(bucket, segments, segment_count, bucket->items[i].start, bucket->items[i].end, &split_id, &merge_id);
                        vg_variate_sv(t_args, split_id, merge_id, &(bucket->items[i]));
                    }
                    break;
                case VG_DIR_OUT:
                    if (bucket->items[i].var == VG_VAR_SNP) {
                        locate_ids(bucket, segments, segment_count, bucket->items[i].start, bucket->items[i].end, &split_id, &merge_id);
                        print_link(split_id, '+', bucket->items[i].id, '+', 0, t_args->out2);
                    } else if (bucket->items[i].var == VG_VAR_INS) {
                        locate_ids(bucket, segments, segment_count, bucket->items[i].start, bucket->items[i].end, &split_id, &merge_id);
                        print_link(split_id, '+', bucket->items[i].id, '+', 0, t_args->out2);
                    } else if (bucket->items[i].var == VG_VAR_INS_SV) {
                        locate_ids(bucket, segments, segment_count, bucket->items[i].start, bucket->items[i].end, &split_id, &merge_id);
                        vg_variate_sv(t_args, split_id, 0, &(bucket->items[i]));
                    } else if (bucket->items[i].var == VG_VAR_ALT) {
                        locate_ids(bucket, segments, segment_count, bucket->items[i].start, bucket->items[i].end, &split_id, &merge_id);
                        print_link(split_id, '+', bucket->items[i].id, '+', 0, t_args->out2);
                    } else if (bucket->items[i].var == VG_VAR_ALT_SV) {
                        locate_ids(bucket, segments, segment_count, bucket->items[i].start, bucket->items[i].end, &split_id, &merge_id);
                        vg_variate_sv(t_args, split_id, 0, &(bucket->items[i]));
                    }
                    break;
                case VG_DIR_INCOMING:
                    locate_ids(bucket, segments, segment_count, bucket->items[i].start, bucket->items[i].end, &split_id, &merge_id);
                    print_link(bucket->items[i].id, '+', merge_id, '+', 0, t_args->out2);
                    break;
                default:
                    fprintf(stderr, "[ERROR] Invalid variation.\n");
                    break;
                }
            }

            // cleanup
            free(bucket->items); free(bucket);
            free(split_points);  free(segments);
        }

        free(batch);
    }

    time_t thread_end;
    time(&thread_end);

    t_args->exec_time = difftime(thread_end, thread_start);
}

// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------
//      MAIN / MASTER
// ------------------------------------------------------------------------------------
// ------------------------------------------------------------------------------------

static inline vg_bucket_batch_t *malloc_vg_bucket_batch(void) {
    vg_bucket_batch_t *batch = (vg_bucket_batch_t *)malloc(sizeof(vg_bucket_batch_t));
    batch->count = 0;
    memset(batch->items, 0, sizeof(batch->items));
    return batch;
}

static inline void reset_vg_bucket_batch(vg_bucket_batch_t *b) {
    b->count = 0;
    memset(b->items, 0, sizeof(b->items));
}

static inline void handle_current_bucket(vg_work_queue_t *queue, vg_bucket_batch_t **batch, vg_core_bucket_t **bucket, vg_queue_sync_t *sync,
                                         struct chr *curr_chr, int chr_idx, int *core_idx, uint64_t *pending_var_ends, int *pending_var_ends_size) {
    // vg_queue_push(queue, *bucket, sync);
    (*batch)->items[(*batch)->count++] = *bucket;

    if ((*batch)->count == VG_BUCKET_BATCH) {
        vg_queue_push(queue, *batch, sync);
        *batch = malloc_vg_bucket_batch();  // start fresh batch
    }

    (*core_idx)++;
    if ((*core_idx) < curr_chr->cores_size) {
        *bucket = malloc_vg_core_bucket(chr_idx, *core_idx, curr_chr->cores[*core_idx].id, curr_chr->cores[(*core_idx) - 1].id);
        add_element_to_bucket(*bucket, pending_var_ends, pending_var_ends_size, curr_chr->cores[*core_idx].start, curr_chr->cores[*core_idx].end);
    } else {
        *bucket = malloc_vg_core_bucket(chr_idx, 0, 0, 0);
    }
}

static inline void flush_batch_if_needed(vg_work_queue_t *queue, vg_bucket_batch_t **batch, vg_queue_sync_t *sync) {
    if (*batch && (*batch)->count > 0) {
        vg_queue_push(queue, *batch, sync);
        *batch = NULL;
    } else if (*batch) {
        free(*batch);
        *batch = NULL;
    }
}

void vg_read_vcf(struct opt_arg *args, struct ref_seq *seqs) {

    printf("[INFO] Processing variations...\n");

    FILE *out_segment = NULL, *out_link = NULL;
    open_files(args, &out_segment, &out_link);

    // create thread arguments
    struct t_arg *t_args = (struct t_arg*)malloc(sizeof(struct t_arg) * args->thread_number);
    pthread_mutex_t queue_mutex   = PTHREAD_MUTEX_INITIALIZER;
    pthread_cond_t cond_not_full  = PTHREAD_COND_INITIALIZER;
    pthread_cond_t cond_not_empty = PTHREAD_COND_INITIALIZER;
    pthread_mutex_init(&queue_mutex, NULL);
    pthread_cond_init(&cond_not_full, NULL);
    pthread_cond_init(&cond_not_empty, NULL);
    int exit_signal = 0;

    vg_queue_sync_t sync = {
        .mutex = &queue_mutex,
        .cond_not_full = &cond_not_full,
        .cond_not_empty = &cond_not_empty,
        .exit_signal = &exit_signal
    };

    vg_work_queue_t queue;
    vg_queue_init(&queue, args->tload_factor * args->thread_number);

    for (int i = 0; i < args->thread_number; i++) {
        char indexed_seg_filename[strlen(args->gfa_path) + 7];
        snprintf(indexed_seg_filename, sizeof(indexed_seg_filename), "%s.s.%d", args->gfa_path, i + 1);
        open_file_w(&(t_args[i].out1), indexed_seg_filename);

        char indexed_lin_filename[strlen(args->gfa_path) + 7];
        snprintf(indexed_lin_filename, sizeof(indexed_lin_filename), "%s.l.%d", args->gfa_path, i + 1);
        open_file_w(&(t_args[i].out2), indexed_lin_filename);

        t_args[i].core_id_index  = ((uint64_t)(i + 1) << 32) + 1;
        t_args[i].thread_id      = i + 1;
        t_args[i].lcp_level      = args->lcp_level;
        t_args[i].is_rgfa        = args->is_rgfa;
        t_args[i].no_overlap     = args->no_overlap;
        t_args[i].seqs           = seqs;
        t_args[i].exec_time      = 0;
        t_args[i].queue          = (void*)&(queue);
        t_args[i].sync           = &sync;
        t_args[i].out_log_mutex  = NULL;
    }

    name_thread("main");

    struct tpool *tm = tpool_create(args->thread_number);

    for (int i = 0; i < args->thread_number; i++) {
        tpool_add_work(tm, vg_read_vcf_thd, t_args + i);
    }

    uint64_t current_size = 1048576;
    char *line = (char *)malloc(current_size);

    FILE *file;
    open_file_r(&file, args->vcf_path);

    int line_count = 0;

    int pending_var_ends_capacity = 256;
    int pending_var_ends_size = 0;
    uint64_t *pending_var_ends = (uint64_t *)malloc(pending_var_ends_capacity * sizeof(uint64_t)); // id+end

    int chr_idx = 0, core_idx = 0, chrom_index = 0;
    struct chr *curr_chr = &(seqs->chrs[chr_idx]);
    vg_bucket_batch_t *batch = malloc_vg_bucket_batch();
    vg_core_bucket_t *bucket = malloc_vg_core_bucket(chr_idx, core_idx, seqs->chrs[chr_idx].cores[core_idx].id, 0);
    curr_chr->ids = (uint64_t **)malloc(curr_chr->cores_size * sizeof(uint64_t *));

    time_t main_start;
    time(&main_start);

    while (fgets(line, current_size, file) != NULL) {
        size_t len = strlen(line);
        int skip_line = 0;
        // read the line and fit it to `line`
        while (len == current_size - 1 && line[len - 1] != '\n') {
            current_size *= 2;
            char *temp_line = (char *)realloc(line, current_size);
            if (!temp_line) {
                fprintf(stderr, "[ERROR] Memory reallocation failed.\n");
                free(line); free(pending_var_ends); free(bucket->items); free(bucket); fclose(file);
                return;
            }
            line = temp_line;

            if (fgets(line + len, current_size - len, file) == NULL) {
                skip_line = 1; break;
            }
            len = strlen(line);
        }
        // validate `line`
        if (skip_line || len < 2 || line[0] == '#') continue;
        if (line[len - 1] == '\n') { line[len - 1] = '\0'; len--; }
        line_count++;

        // parse the `line`
        char *chrom, *index, *id, *ref, *alt;
        size_t offset;

        char *saveptr;
        chrom = strtok_r(line, "\t", &saveptr); // get chromosome name
        if (strcmp(chrom, seqs->chrs[chrom_index].seq_name) != 0) {
            if (chrom_index + 1 < seqs->size && strcmp(chrom, seqs->chrs[chrom_index + 1].seq_name) == 0) { // if it continues with right next chromosome
                chrom_index++;
            } else { // search right chromosome from the beggining
                chrom_index = -1;
                for (int i = chr_idx; i < seqs->size; i++) {
                    if (strcmp(chrom, seqs->chrs[i].seq_name) == 0) { chrom_index = i; break; }
                }
                // if the chromosome is not found (rare or if VCF is not compatible with fasta)
                if (chrom_index == -1) continue;
            }
        }
        
        index   = strtok_r(NULL, "\t", &saveptr);   // get index  
        offset  = strtol(index, NULL, 10) - 1;      // get offset
        id      = strtok_r(NULL, "\t", &saveptr);   // get ID
        ref     = strtok_r(NULL, "\t", &saveptr);   // get REF
        alt     = strtok_r(NULL, "\t", &saveptr);   // get ALT alleles

        // If we move to next lcp core, push array if there are elements and create new array
        if (chrom_index == chr_idx && curr_chr->cores[core_idx].end <= offset) {
            while (core_idx < curr_chr->cores_size && curr_chr->cores[core_idx].end <= offset) {
                handle_current_bucket(&queue, &batch, &bucket, &sync, curr_chr, chr_idx, &core_idx, pending_var_ends, &pending_var_ends_size);
            }
        } else if (chrom_index != chr_idx) {
            // it seems that the vcf file moved to new chromosome. then, print remaining lcp cores on prev chrom
            while (core_idx < curr_chr->cores_size) {
                handle_current_bucket(&queue, &batch, &bucket, &sync, curr_chr, chr_idx, &core_idx, pending_var_ends, &pending_var_ends_size);
            }
            chr_idx++;
            
            // if there is a chromosomal jump (e.g., from chr1 to chr4), print chr2 and chr3
            while (chr_idx < chrom_index) {
                vg_print_seq(&(seqs->chrs[chr_idx]), args->is_rgfa, out_segment, out_link);
                chr_idx++;
            }
            
            // reset chromosome and index info as it is a new chromosome
            chr_idx = chrom_index;
            core_idx = 0;
            curr_chr = &(seqs->chrs[chr_idx]);
            curr_chr->ids = (uint64_t **)malloc(curr_chr->cores_size * sizeof(uint64_t *));
            
            // move bucket data to correct position
            while (core_idx < curr_chr->cores_size && curr_chr->cores[core_idx].end <= offset) {
                vg_print_core_as_is(curr_chr, chr_idx, core_idx, seqs, args->is_rgfa, out_segment, out_link);
                core_idx++;
            }

            // reset bucket data
            bucket->chr_idx  = chrom_index;
            bucket->core_idx = core_idx;
            bucket->prev_id  = core_idx ? curr_chr->cores[core_idx - 1].id : 0;
            bucket->curr_id  = curr_chr->cores[core_idx].id;
        }

        // ALT can be multi-allelic; store one element per ALT if you want
        size_t rlen = strlen(ref);
        size_t alen = strlen(alt);
        int order = 0;

        if (rlen > 1 && alen == 1) {
            char *ref_saveptr;
            char *ref_token = strtok_r(ref, ",", &ref_saveptr); // split REF alleles by comma (it is rare but in case it happens)
            while (ref_token != NULL) {
                if (!check_vg_items(bucket)) break;

                size_t tlen = strlen(ref_token);
                
                if (offset + tlen < curr_chr->cores[core_idx].end) { // DEL inside
                    bucket->items[bucket->size] = (vg_element_t){VG_DIR_IN, VG_VAR_DEL, 0, offset + 1, offset + tlen, NULL, NULL, order};
                    bucket->size++;
                } else if (offset + 1 < curr_chr->cores[core_idx].end) { // it starts inside
                    add_pending_var_end(&pending_var_ends, &pending_var_ends_size, &pending_var_ends_capacity, args->core_id_index, offset + tlen);
                    bucket->items[bucket->size] = (vg_element_t){VG_DIR_OUT, VG_VAR_DEL, args->core_id_index, offset + 1, 0xFFFFFFFFFFFFFFFF, NULL, NULL, order};
                    bucket->size++;
                    args->core_id_index++;
                } else {
                    add_pending_var_end(&pending_var_ends, &pending_var_ends_size, &pending_var_ends_capacity, curr_chr->cores[core_idx].id, offset + tlen);
                }
                
                ref_token = strtok_r(NULL, ",", &ref_saveptr);
                order++;
            }
            continue;
        }

        char *alt_saveptr;
        char *alt_token = strtok_r(alt, ",", &alt_saveptr); // split ALT alleles by comma
        while (alt_token != NULL) {
            if (!check_vg_items(bucket)) break;

            size_t tlen = strlen(alt_token);

            if (rlen == 1 && tlen == 1) { // SNP
                if (offset + 1 < curr_chr->cores[core_idx].end) {
                    print_seq_vg(args->core_id_index, alt_token, 1, id, order, offset, 1, args->is_rgfa, out_segment);
                    bucket->items[bucket->size] = (vg_element_t){VG_DIR_IN, VG_VAR_SNP, args->core_id_index, offset, offset + 1, NULL, NULL, order}; // id assigned for segment
                } else {
                    add_pending_var_end(&pending_var_ends, &pending_var_ends_size, &pending_var_ends_capacity, args->core_id_index, offset + 1);
                    print_seq_vg(args->core_id_index, alt_token, 1, id, order, offset, 1, args->is_rgfa, out_segment);
                    bucket->items[bucket->size] = (vg_element_t){VG_DIR_OUT, VG_VAR_SNP, args->core_id_index, offset, 0xFFFFFFFFFFFFFFFF, NULL, NULL, order}; // id assigned for segment
                }
                bucket->size++;
                args->core_id_index++;
            } else if (1 == rlen) { // INS
                // Small insertion
                if (tlen / 2 < curr_chr->cores[core_idx].end - curr_chr->cores[core_idx].start) {
                    print_seq_vg(args->core_id_index, alt_token + 1, tlen - 1, id, order, offset, 1, args->is_rgfa, out_segment);
                    if (offset + 1 < curr_chr->cores[core_idx].end) {
                        bucket->items[bucket->size] = (vg_element_t){VG_DIR_IN, VG_VAR_INS, args->core_id_index, offset + 1, offset + 1, NULL, NULL, order}; // id assigned for segment         
                    } else {
                        add_pending_var_end(&pending_var_ends, &pending_var_ends_size, &pending_var_ends_capacity, args->core_id_index, offset + 1);
                        bucket->items[bucket->size] = (vg_element_t){VG_DIR_OUT, VG_VAR_INS, args->core_id_index, offset + 1, 0xFFFFFFFFFFFFFFFF, NULL, NULL, order}; // id assigned for segment
                    }
                    bucket->size++;
                    args->core_id_index++;
                } else {  // Large INS, to be processed with LCP
                    char *alt_token_copy = strdup(alt_token);
                    char *seq_id = strdup(id);
                    if (offset + 1 < curr_chr->cores[core_idx].end) { // if inside of the lcp core
                        bucket->items[bucket->size] = (vg_element_t){VG_DIR_IN, VG_VAR_INS_SV, args->core_id_index, offset + 1, offset + 1, alt_token_copy, seq_id, order};
                    } else { // if in the edge of the end of the lcp core
                        add_pending_var_end(&pending_var_ends, &pending_var_ends_size, &pending_var_ends_capacity, args->core_id_index, offset + 1);
                        bucket->items[bucket->size] = (vg_element_t){VG_DIR_OUT, VG_VAR_INS_SV, args->core_id_index, offset + 1, 0xFFFFFFFFFFFFFFFF, alt_token_copy, seq_id, order};
                    }
                    bucket->size++;
                    args->core_id_index++;
                }
            } else { // ALT
                if (tlen / 2 < curr_chr->cores[core_idx].end - curr_chr->cores[core_idx].start) { // alteration, simply print the underling string
                    print_seq_vg(args->core_id_index, alt_token, tlen, id, order, offset, 1, args->is_rgfa, out_segment);
                    if (offset + rlen < curr_chr->cores[core_idx].end) {
                        bucket->items[bucket->size] = (vg_element_t){VG_DIR_IN, VG_VAR_ALT, args->core_id_index, offset, offset + rlen, NULL, NULL, order};
                    } else {
                        add_pending_var_end(&pending_var_ends, &pending_var_ends_size, &pending_var_ends_capacity, args->core_id_index, offset + rlen);
                        bucket->items[bucket->size] = (vg_element_t){VG_DIR_OUT, VG_VAR_ALT, args->core_id_index, offset, 0xFFFFFFFFFFFFFFFF, NULL, NULL, order};
                    }
                    bucket->size++;
                    args->core_id_index++;
                } else { // check if it the alt_token requires LCP processing
                    char *alt_token_copy = strdup(alt_token);
                    char *seq_id = strdup(id);
                    if (offset + rlen < curr_chr->cores[core_idx].end) {
                        bucket->items[bucket->size] = (vg_element_t){VG_DIR_IN, VG_VAR_ALT_SV, args->core_id_index, offset, offset + rlen, alt_token_copy, seq_id, order};
                    } else {
                        add_pending_var_end(&pending_var_ends, &pending_var_ends_size, &pending_var_ends_capacity, args->core_id_index, offset + rlen);
                        bucket->items[bucket->size] = (vg_element_t){VG_DIR_OUT, VG_VAR_ALT_SV, args->core_id_index, offset, 0xFFFFFFFFFFFFFFFF, alt_token_copy, seq_id, order};
                    }
                    bucket->size++;
                    args->core_id_index++;
                }
            }
            alt_token = strtok_r(NULL, ",", &alt_saveptr);
            order++;
        }
    }

    // handle remaining LCP cores
    while (core_idx < curr_chr->cores_size) {
        // if there is anything to push as a job into the pool, then push it (previous core's data)
        handle_current_bucket(&queue, &batch, &bucket, &sync, curr_chr, chr_idx, &core_idx, pending_var_ends, &pending_var_ends_size);
    }
    
    if (batch->count) {
        flush_batch_if_needed(&queue, &batch, &sync);
    } else {
        free(batch);
    }

    chr_idx++;
    // print remaining chromosomes if any
    while (chr_idx < seqs->size) {
        vg_print_seq(&(seqs->chrs[chr_idx]), args->is_rgfa, out_segment, out_link);
        chr_idx++;
    }
    fclose(file);

    pthread_mutex_lock(sync.mutex);
    while (queue.size > 0) { 
        pthread_cond_wait(sync.cond_not_full, sync.mutex);
    }
    pthread_mutex_unlock(sync.mutex);

    exit_signal = 1;
    pthread_cond_broadcast(sync.cond_not_empty);
    tpool_wait(tm);
    tpool_destroy(tm);
    pthread_mutex_destroy(sync.mutex);
    pthread_cond_destroy(sync.cond_not_full);
    pthread_cond_destroy(sync.cond_not_empty);

    time_t main_end;
    time(&main_end);

    for (int i = 0; i < args->thread_number; i++) {
        fclose(t_args[i].out1);
        fclose(t_args[i].out2);
    }
    free(t_args);

    printf("[INFO] VCF processing completed in %0.2f sec.\n", difftime(main_end, main_start));

    while (queue.size) {
        fprintf(stderr, "[WARN] Left work in the queue.\n"); // should not happen
        vg_bucket_batch_t *batch = queue.batch[queue.front];
        queue.front = (queue.front + 1) % queue.capacity;
        queue.size--;
        if (batch) {
            for (int j = 0; j < batch->count; j++) {
                vg_core_bucket_t *bucket = batch->items[j];
                for (int i = 0; i < bucket->size; i++) {
                    if (bucket->items[i].seq != NULL)    free(bucket->items[i].seq);
                    if (bucket->items[i].seq_id != NULL) free(bucket->items[i].seq_id);
                }
                free(bucket->items);
                free(bucket);
            }
            free(batch);
        }
    }
    free(queue.batch);
    free(line);

    free(pending_var_ends);
    
    // print path
    FILE *out_path;
    char path_filename[strlen(args->gfa_path)+5];
    snprintf(path_filename, sizeof(path_filename), "%s.p", args->gfa_path);
    out_path = fopen(path_filename, "w");
    if (out_path == NULL) {
        fprintf(stderr, "[ERROR] Couldn't open path file\n");
        exit(EXIT_FAILURE);
    }

    print_path(seqs, out_path);
    fclose(out_path);
    fclose(out_segment);
    fclose(out_link);
}