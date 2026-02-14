#ifndef __STRUCT_DEF_H__
#define __STRUCT_DEF_H__

#include <stdint.h>
#include <stdio.h>
#include <pthread.h>

#define THREAD_POOL_FACTOR 2

typedef enum {
    VG,
    VGX,
    LDBG
} program_mode;

struct opt_arg {
	char *fasta_path;		/** Path to the input FASTA file. */
	char *fasta_fai_path;	/** Path to the input FASTA  index file. */
	char *vcf_path;			/** Path to the input VCF file. */
	char *gfa_path; 		/** Path to the output rGFA/GFA file. */
    char *prefix;           /** Prefix to the files */
    program_mode program;   /** Program mode. */
	uint64_t core_id_index; /** Global id index for LCP cores. */
	int lcp_level;			/** The LCP level to be used. */
    int thread_number;      /** The thread number. */
    int failed_var_count;   /** Total number of lines failed to variate. */
    int invalid_line_count; /** Total number of lines that are invalid in VCF file. */
	int bubble_count;	 	/** Number of bubbles created in the graph. */
	int is_rgfa;			/** Boolean argument to output rGFA or GFA. */
    int no_overlap;         /** Boolean argument to decide whether allow overlap. */
    int skip_masked;        /** Boolean argument to decide whether include invalid chars (N) to the output. */
    int tload_factor;       /** Thread pool element storage capacity factor to the tread number. */
    int verbose;            /** Verbose. */
};

struct simple_core {
	uint64_t id;    /** Core id given by lcpan.*/
	uint64_t start; /** Start index of core. */
	uint64_t end;   /** End index of core. */
};

struct chr {
    char *seq_name;            /** Chromosome name */
    uint64_t global_index;     /** Global Start index (cumulative index from previous chrs). */
    int seq_size;              /** Chromosome size */
    char *seq;                 /** Chromosome Sequence */
	int cores_size;			   /** LCP cores count in cores arrat */
	struct simple_core *cores; /** LCP (ordered) cores in the chromosome */ 
    uint64_t **ids;            /** IDs of sub-segments splitted in the segment (needed for vg-path). */
};

struct ref_seq {
	int size;         /** Number of chromosomes. */
	struct chr *chrs; /** Array of chromosomes */
};

struct line_queue {
    char **lines; /** Queue to store lines extracted from VCF file for threads. */
    int size;     /** The size of the queue. */
    int capacity; /** Capacity of the queue. */
    int front;    /** The index for the pushing point. */
    int rear;     /** The index for the popping point. */
};

typedef enum {
    VG_DIR_IN,
    VG_DIR_OUT,
    VG_DIR_INCOMING
} vg_direction_t;

typedef enum {
    VG_VAR_SNP,
    VG_VAR_INS,
    VG_VAR_INS_SV,
    VG_VAR_DEL,
    VG_VAR_ALT,
    VG_VAR_ALT_SV,
    VG_VAR_NONE
} vg_variant_t;

typedef struct {
    vg_direction_t dir;     /** variation direction */
    vg_variant_t   var;     /** variation type */
    uint64_t id;            /** id of the lcp core */
    uint64_t start;         /** start of the lcp core */
    uint64_t end;           /** end of the lcp core */
    char *seq;              /** the chromosome seqeunce that lcp core lies */
    char *seq_id;           /** the name of the chromosome sequence */
    int order;              /** the variation's index. there might be multiple haplotides in the same index. helps to determine the order */
} vg_element_t;

typedef struct {
    int chr_idx;                  /** Chromosome index. */
    int core_idx;                 /** LCP core index in chromosome. */
    
    int capacity;                 /** Capacity of variation array. */
    int size;                     /** Size of variation array. */
    
    uint64_t prev_id;             /** previous segment's id. */
    uint64_t curr_id;             /** current segment's id. */
    
    vg_element_t *items;          /** Pointer to variations array. */
} vg_core_bucket_t;

typedef struct {
    vg_core_bucket_t **items;   /** Queue to data of the variations to be processed. */
    int size;                   /** The size of the queue. */
    int capacity;               /** Capacity of the queue. */
    int front;                  /** The index for the pushing point. */
    int rear;                   /** The index for the popping point. */
} vg_work_queue_t;

struct t_arg {
    uint64_t core_id_index;
    int thread_id;
    int lcp_level;
	int is_rgfa;
    int no_overlap;
    struct ref_seq *seqs;
    int failed_var_count;
    int invalid_line_count;
    int bubble_count;
    double exec_time;
    FILE *out1;
    FILE *out2;
    void *queue;
    pthread_mutex_t *queue_mutex;
    pthread_mutex_t *out_log_mutex;
    pthread_cond_t *cond_not_full;
    pthread_cond_t *cond_not_empty;
    int *exit_signal;
};

#endif