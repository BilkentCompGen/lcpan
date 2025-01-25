#ifndef __STRUCT_DEF_H__
#define __STRUCT_DEF_H__

#include <stdint.h>

struct opt_arg {
	char *fasta_path;		/**< Path to the input FASTA file. */
	char *fasta_fai_path;	/**< Path to the input FASTA  index file. */
	char *vcf_path;			/**< Path to the input VCF file. */
	char *gfa_path; 		/**< Path to the output rGFA/GFA file. */
	uint64_t core_id_index; /**< Global id index for LCP cores. */
	int lcp_level;			/**< The LCP level to be used. */
	int bubble_count;	 	/**< Number of bubbles created in the graph. */
	int is_rgfa;			/**< Boolean argument to output rGFA or GFA. */
    int no_overlap;         /**< Boolean argument to decide whether allow overlap. */
};

struct simple_core {
	uint64_t id;
	uint64_t start;
	uint64_t end;
	uint32_t label;
};

struct chr {
    char *seq_name;             /**< Chromosome name */
    int seq_size;               /**< Chromosome size */
    char *seq;                  /**< Chromosome Sequence */
	int cores_size;			    /**< LCP cores count in cores arrat */ 
	struct simple_core *cores;	/**< LCP (ordered) cores in the chromosome */ 
};

struct ref_seq {
	int size;
	struct chr *chrs;
};

#endif