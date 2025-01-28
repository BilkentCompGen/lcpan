#ifndef __STRUCT_DEF_H__
#define __STRUCT_DEF_H__

#include <stdint.h>
#include <stdio.h>
#include <pthread.h>

struct opt_arg {
	char *fasta_path;		/**< Path to the input FASTA file. */
	char *fasta_fai_path;	/**< Path to the input FASTA  index file. */
	char *vcf_path;			/**< Path to the input VCF file. */
	char *gfa_path; 		/**< Path to the output rGFA/GFA file. */
	uint64_t core_id_index; /**< Global id index for LCP cores. */
	int lcp_level;			/**< The LCP level to be used. */
    int thread_number;      /**< The thread number. */
	int bubble_count;	 	/**< Number of bubbles created in the graph. */
	int is_rgfa;			/**< Boolean argument to output rGFA or GFA. */
    int no_overlap;         /**< Boolean argument to decide whether allow overlap. */
};

struct simple_core {
	uint64_t id;    /**< Core id given by lcpan.*/
	uint64_t start; /**< Start index of core. */
	uint64_t end;   /**< End index of core. */
	uint32_t label; /**< Core id given by lcptools. */
};

struct chr {
    char *seq_name;             /**< Chromosome name */
    int seq_size;               /**< Chromosome size */
    char *seq;                  /**< Chromosome Sequence */
	int cores_size;			    /**< LCP cores count in cores arrat */ 
	struct simple_core *cores;	/**< LCP (ordered) cores in the chromosome */ 
};

struct ref_seq {
	int size;           /**< Number of chromosomes. */
	struct chr *chrs;   /**< Array of chromosomes */
};

struct line_queue {
    char **lines;   /**< Queue to store lines extracted from VCF file for threads. */
    int size;       /**< The size of the queue. */
    int capacity;   /**< Capacity of the queue. */
    int front;      /**< The index for the pushing point. */
    int rear;       /**< The index for the popping point. */
};

struct t_arg {
    struct opt_arg *args;
    struct ref_seq *seqs;
    int failed_var_count;
    int invalid_line_count;
    FILE *out;
    FILE *out_err;
    struct line_queue *queue;
    pthread_mutex_t *mutex;
    pthread_mutex_t *out_err_mutex;
    pthread_mutex_t *core_id_mutex;
    pthread_cond_t *cond_not_full;
    pthread_cond_t *cond_not_empty;
    int *exit_signal;
};

#endif