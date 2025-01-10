#ifndef __STRUCT_DEF_H__
#define __STRUCT_DEF_H__

#include <stdint.h>

struct opt_arg {
	char *fasta_path;		/**< Path to the input FASTA file. */
	char *fasta_fai_path;	/**< Path to the input FASTA  index file. */
	char *vcf_path;			/**< Path to the input VCF file. */
	char *rgfa_path; 		/**< Path to the output rGFA file. */
	uint64_t core_id_index; /**< Global id index for LCP cores. */
	int lcp_level;			/**< The LCP level to be used. */
	int bubble_count;	 	/**< Number of bubbles created in the graph. */
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

// std::vector<lcp::lps*> variated_strs;

// // Sequence graph node structure
// /**
//  * @struct core_node
//  * Represents a node in the sequence graph.
//  */
// struct core_node {
//     std::vector<core_node*> next; /**< Pointers to adjacent nodes. */
//     lcp::core* core_value = nullptr;
//     std::string id;                /**< Unique identifier for the node. */
//     bool end_flag = false;       /**< Indicates if the node is the end of the sequence. */
//     std::vector<int> SN_ids;
//     int SO;
//     bool rank;
// };

// // Data structures for the rGFA format
// /**
//  * @struct segment
//  * Represents a segment in the rGFA format.
//  */
// struct segment {
//     core_node* core;       	/**< The nucleotide sequence of the segment. */
//     std::string seg_name;   /**< Identifier for the segment. */
// };

// /**
//  * @struct g_link
//  * Represents a graph link in the rGFA format.
//  */
// struct g_link {
//     std::string source;         /**< Source segment ID. */
//     bool source_orient;         /**< Orientation of the source segment. */
//     std::string target;         /**< Target segment ID. */
//     bool target_orient;         /**< Orientation of the target segment. */
//     u_int cigar;                /**< Edit distance information in CIGAR format. */
// };

// std::vector<segment*> segments;
// std::vector<g_link*> links;
// std::string ref_path = "P ref ";

// int curr_indx = 0;
// std::vector<core_node*> indices;

// // Sequence graph class
// /**
//  * @class sequence_graph
//  * Represents the sequence graph for the genome.
//  */
// class sequence_graph {
//     public:
//         core_node* head = nullptr; /**< Pointer to the first node in the graph. */
//         int size = 0;              /**< Total number of nodes in the graph. */
// };


#endif