#define STATS

/**
 * @file vg_lcp.cpp
 * @brief Implementation of a tool for constructing sequence graphs from genomic data.
 *
 * This file reads FASTA input and identifies its LCP (Locally Consistent Parsing) cores. 
 * It processes variations from a VCF file by locating the affected regions on the 
 * original sequence, identifying the LCP cores of those regions, aligning the original 
 * and varied cores, and creating bubbles to represent the differences. 
 * After applying all variations, it constructs a variation graph and outputs it in 
 * rGFA format, representing genomic sequences and their relationships for graph-based 
 * genome assembly and analysis.
 * @note The `STATS` macro enables performance analytics by tracking core start and 
 *       end indices.
 *
 * Structures used:
 * - `segment`: Represents a segment in the rGFA format.
 * - `g_link`: Represents a graph link in the rGFA format.
 * - `core_node`: Represents a node in the sequence graph with pointers to other nodes 
 *   and its LCP core value.
 *
 * Utility functions:
 * - `core_to_seq()`: Converts a core string representation to a nucleotide sequence.
 * - `dfs()`: Performs a depth-first traversal of the sequence graph to generate graph links.
 * - `find_boundaries()`: Determines the boundaries of a variation in the sequence graph.
 * - `align_variation()`: Aligns a variation against the sequence graph to identify differences.
 * - `construct_variated_seq()`: Constructs a modified sequence incorporating a variation.
 * - `variate()`: Applies variations to the sequence graph, creating bubble structures.
 *
 * Command-line options:
 * - `-f`: Specifies the input FASTA file path.
 * - `-v`: Specifies the input VCF file path.
 * - `-r`: Specifies the output rGFA file path.
 * - `-l`: Specifies the LCP level for parsing.
 */


#include <string.h>
#include <sstream>
#include <iostream>
#include <vector>
#include <unordered_set>
#include <fstream>
#include <cmath>
#include <unistd.h>

#include "lps.h"
#include "fasta_reader.h"
#include "vcf_parser.h"
#include <limits>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

// Global variables and constants
int global_no_of_bubbles = 0; /**< Number of bubbles created in the graph. */
std::string fasta_path;       /**< Path to the input FASTA file. */
std::string vcf_path;         /**< Path to the input VCF file. */
std::string rgfa_path;        /**< Path to the output rGFA file. */
int level;                    /**< LCP level for parsing. */
int offsets[] = {1, 4, 11, 27}; /**< Offsets for LCP parsing levels obtained experimentally. */
std::vector<std::string> chrmsms;
int idx_size = 100;

std::vector<lcp::lps*> variated_strs;

// Sequence graph node structure
/**
 * @struct core_node
 * Represents a node in the sequence graph.
 */
typedef struct core_node {
    std::vector<core_node*> next; /**< Pointers to adjacent nodes. */
    lcp::core* core_value = nullptr;
    std::string id;                /**< Unique identifier for the node. */
    bool end_flag = false;       /**< Indicates if the node is the end of the sequence. */
    std::vector<int> SN_ids;
    int SO;
    bool rank;
} core_node;

// Data structures for the rGFA format
/**
 * @struct segment
 * Represents a segment in the rGFA format.
 */
typedef struct segment {
    core_node* core;       /**< The nucleotide sequence of the segment. */
    std::string seg_name;   /**< Identifier for the segment. */
} segment;

/**
 * @struct g_link
 * Represents a graph link in the rGFA format.
 */
typedef struct g_link {
    std::string source;         /**< Source segment ID. */
    bool source_orient;         /**< Orientation of the source segment. */
    std::string target;         /**< Target segment ID. */
    bool target_orient;         /**< Orientation of the target segment. */
    u_int cigar;                /**< Edit distance information in CIGAR format. */
} g_link;

std::vector<segment*> segments;
std::vector<g_link*> links;
std::string ref_path = "P ref ";

int curr_indx = 0;
std::vector<core_node*> indices;

// Sequence graph class
/**
 * @class sequence_graph
 * Represents the sequence graph for the genome.
 */
class sequence_graph {
    public:
        core_node* head = nullptr; /**< Pointer to the first node in the graph. */
        int size = 0;              /**< Total number of nodes in the graph. */
};

// Function implementations
/**
 * @brief Converts a core string representation to a nucleotide sequence.
 * @param core The core string to be converted.
 * @return The corresponding nucleotide sequence.
 */
std::string core_to_seq(std::string core) {
    std::string seq = "";
    for (int i = 0; i < (int) (core.size() - 1); i += 2) {
        std::string char_rep = core.substr(i, 2);
        if (char_rep.compare("00") == 0) {
            seq += "A";
        } else if (char_rep.compare("01") == 0) {
            seq += "C";
        } else if (char_rep.compare("10") == 0) {
            seq += "G";
        } else {
            seq += "T";
        }
    }
    return seq;
}

/**
 * @brief Performs a depth-first traversal of the sequence graph to securely delete the core_node objects.
 * @param node The starting node for traversal.
 * @param visited A set to track visited nodes.
**/
void delete_core_nodes(core_node* node, std::unordered_set<core_node*>& visited) {
    if (node == nullptr || visited.count(node)) {
        return;
    }

    visited.insert(node);

    // Recursively delete all connected nodes
    for (core_node* next_node : node->next) {
        delete_core_nodes(next_node, visited);
    }

    // Delete the current node after visiting its neighbors
    delete node;
}


/**
 * @brief Calculates the offset for LCP parsing based on the sequence length.
 * @param n The sequence length.
 * @return The calculated offset.
 */
int calculate_offset(int n) {
    return ceil(3.6 * pow((1 / 0.45), (n-1)));
}

/**
 * @brief Performs a depth-first traversal of the sequence graph.
 * @param node The starting node for traversal.
 * @param visited A set to track visited nodes.
 * @param in_bubble Indicates if the traversal is within a bubble structure.
 */
void dfs(core_node* node, std::unordered_set<core_node*>& visited) {
    if (node == nullptr || visited.count(node)) {
        return;
    }

    visited.insert(node);

    std::ostringstream oss;
    segment* seg = new segment();
    seg->core = node;
    seg->seg_name = node->id;

    segments.push_back(seg);

    if (node->id.at(0) == 'o') ref_path += (node->id + "+, ");

    for (int i = 0; i < (int) (node->next.size()); i++) {

        g_link* lnk = new g_link();
        lnk->source = node->id;
        if (node->next.at(i)->id.at(0) == 'b' && node->next.at(i)->id.at(1) == 'h') {
            lnk->target = "00";
        } else {
            lnk->target = node->next.at(i)->id;
        }
        lnk->source_orient = 1;
        lnk->target_orient = 1;
        lnk->cigar = 0;
        links.push_back(lnk);


        core_node* next_node = node->next.at(i);
        if (node->id.at(0) == 'b' && node->next.at(i)->id.at(0) == 'o') {
            std::cout << node->core_value << ":" << node->id << ":";
            for (int i = 0; i < (int) (node->next.size()); i++) {
                std::cout << node->next.at(i)->id << ", ";
            }
            std::cout << "\n";
            continue;
        } else if (node->id.at(0) == 'o' && node->next.at(i)->id.at(0) == 'b') {
            std::cout << node->core_value << ":" << node->id << ":";
            for (int i = 0; i < (int) (node->next.size()); i++) {
                std::cout << node->next.at(i)->id << ", ";
            }
            std::cout << "\n";
            dfs(next_node, visited);
        } else if (node->id.at(0) == 'o' && node->next.at(i)->id.at(0) == 'o') {
            dfs(next_node, visited);
        } else {
            std::cout << node->core_value << ":" << node->id << ":";
            for (int i = 0; i < (int) (node->next.size()); i++) {
                std::cout << node->next.at(i)->id << ", ";
            }
            std::cout << "\n";
            dfs(next_node, visited);
        }
    }
}

/**
 * @brief Prints the sequence graph and writes it in rGFA format.
 * @param graph The sequence graph to be printed.
 */
void print_seq(sequence_graph& graph, std::string rgfa_path, std::string c) {
    std::unordered_set<core_node*> visited;

    std::cout << "\n********************************************\n" << std::endl;

    if (graph.head != nullptr) {
        dfs(graph.head, visited);
    }

    std::cout << "\n\n********************************************\n" << std::endl;

    std::ofstream out_file(rgfa_path);
    for (segment* s : segments) {
        out_file << "S\t" << s->seg_name << "\t" << s->core->core_value << "\tSN:Z:";
        
        for (int i : s->core->SN_ids) {
            if (i == -1) out_file << c << ",";
            else {
                out_file << chrmsms.at(i) <<",";
            }
        }

            out_file << "\tSO:i:" << s->core->SO << "\tSR:i:" << (s->core->rank ? 0 : 1);
        out_file << std::endl;
    }

    for (g_link* l : links) {
        if (l->target.at(0) == '0' && l->target.at(1) == '0') continue;
        out_file << "L\t" << l->source << "\t" << (l->source_orient ? "+\t" : "-\t") <<
            l->target << "\t" << (l->target_orient ? "+\t" : "-\t") << l->cigar << "M";
        out_file << std::endl;
    }

    out_file << ref_path;
    out_file.close();
}

/**
 * @brief Determines the boundaries of a variation in the sequence graph.
 * @param start_loc The start position of the variation.
 * @param end_loc The end position of the variation.
 * @return Array containing the start index, end index, start core number, and end core number.
 */
int* find_boundaries(int start_loc, int end_loc) {

    int start_index; // indices in the original sequence
    int end_index; 
    int start_core_number = 0; // start and end core numbers fixing the first possible changed node as 0
    int end_core_number = 0;

    int offset = level <= 4 ? (int)(offsets[level - 1] * 1.5) : calculate_offset(level) * 1.5;

    start_loc = MAX(0, start_loc - offset);
    end_loc = MIN(end_loc + offset, std::numeric_limits<int>::max());

    // traverse to find the starting core
    core_node* curr_start;
    core_node* prev_start;

    int indices_idx = start_loc / idx_size;
    if (indices_idx >= (int) indices.size()) indices_idx = indices.size() - 1;
    core_node* temp = indices.at(indices_idx);
    curr_start = ((int) temp->core_value->start <= start_loc || indices_idx == 0) ? temp : indices.at(indices_idx - 1);
    prev_start = curr_start;

    bool end_case = false;

    while ((int) (curr_start->core_value->end) < start_loc) {
        if (curr_start->end_flag) {
            end_case = true;
            break;
        }
        prev_start = curr_start;
        curr_start = curr_start->next.at(0);
        start_core_number++;
    }

    // traverse to find the ending core
    core_node* curr_end = curr_start;
    end_core_number = start_core_number - 1;

    while ((int) (curr_end->core_value->start) <= end_loc && !end_case) {
        if (curr_end->end_flag) {
            end_case = true;
            break;
        }
        curr_end = curr_end->next.at(0);
        end_core_number++;
    }

    // set the start and the end indices in the original sequence
    start_index = prev_start->core_value->start; 
    end_index = curr_end->core_value->end;

    int* boundaries = new int[4]{start_index, end_index, start_core_number, end_core_number};
    return boundaries;
}

/**
 * @brief Aligns a variation with the original sequence graph to determine differences.
 * @param starting_core The starting core node in the original sequence graph.
 * @param var_str The variation sequence represented as an LCP structure.
 * @param strt_core_in_org The index of the starting core in the original sequence.
 * @param end_core_in_org The index of the ending core in the original sequence.
 * @param first_difference Outputs the index of the first differing core.
 * @param last_difference_org Outputs the last differing core in the original sequence.
 * @param last_difference_var Outputs the last differing core in the variation sequence.
 * @return The core node where alignment starts.
 */
core_node* align_variation(core_node* starting_core, lcp::lps* var_str,  
    int strt_core_in_org, int end_core_in_org,
    int& first_difference, int& last_difference_org, int& last_difference_var) {
    

    if (var_str->size() <= 0) {
	first_difference = -1;
	last_difference_org = -1;
	last_difference_var = -1;
	return starting_core;
	
    }	
    // find when the two cores set are different for the first time
    first_difference = 0;
    core_node* curr = starting_core->next.at(0);
    core_node* prev_start = starting_core;
    bool diff = true;

    curr = prev_start->next.at(0);

    while (curr != nullptr && !curr->end_flag && first_difference < (int) (var_str->size() - 1) && 
        curr->core_value->label == var_str->cores->at(first_difference).label) {
            
       
        diff = false;

        first_difference++;
        prev_start = curr;
        if (curr->next.size() > 0)
            curr = curr->next.at(0);
        else 
            curr = nullptr;
    }


    // create stacks with enough spaces for the original and the varaited sequences
    core_node* original_stack [end_core_in_org - strt_core_in_org + 1];
    core_node* variated_stack [var_str->size()];
    int original_stack_size = 0;
    int variated_stack_size = 0;

    core_node* return_core = prev_start;
    curr = starting_core->next.at(0);

    // add the cores to the stacks
    for (int i = strt_core_in_org; i <= end_core_in_org; i++) {
        core_node* new_cn = new core_node();
        new_cn->core_value = curr->core_value;
        original_stack[original_stack_size++] = new_cn;
        if (curr->next.size() > 0)
            curr = curr->next.at(0);
    }

    for (int i = 0; i < (int) (var_str->cores->size()); i++) {
        core_node* new_cn = new core_node();
        new_cn->core_value = &var_str->cores->at(i);
        variated_stack[variated_stack_size++] = new_cn;
    }

    // find the first node they are different from the end
    last_difference_org = end_core_in_org - strt_core_in_org;
    last_difference_var = var_str->cores->size()-1;
    
            
    while (original_stack_size > 0 && variated_stack_size > 0 && 
        original_stack[--original_stack_size]->core_value->label == variated_stack[--variated_stack_size]->core_value->label) {
            last_difference_org--;
            last_difference_var--;
    }

    if (diff) first_difference = -1;

    for (core_node* o : original_stack) {
        delete o;
    }

    for (core_node* v : variated_stack) {
        delete v;
    }

    return return_core;
}


/**
 * @brief Constructs the sequence with the applied variation.
 * @param sequence The original sequence.
 * @param variation The variation to apply (insert, delete, or substitute).
 * @param start_loc The starting position for the variation.
 * @param end_loc The ending position for the variation.
 * @param variation_type The type of variation (0: insertion, 1: substitution, 2: deletion).
 * @return The newly constructed variated sequence.
 */
std::string construct_variated_seq(std::string &sequence, std::string variation, 
    int start_loc, int end_loc, int variation_type) {

    int* boundaries = find_boundaries(start_loc, end_loc);
    
    // construct the new varaited string and find its cores
    std::string before, after;
    std::string variated_seq;

    switch (variation_type)
    {
    case 0: //insertion
        before = sequence.substr(boundaries[0] - 1, start_loc - boundaries[0]);
        after = sequence.substr(end_loc - 1, boundaries[1] - end_loc + 1);
        variated_seq = before + variation + after;
        break;
    case 1: //diff
        
        before = sequence.substr(boundaries[0] - 1, start_loc - boundaries[0]);
        after = sequence.substr(end_loc, boundaries[1] - end_loc);
        variated_seq = before + variation + after;
        break;
    case 2: //diff
        // todo: change length independent from level

        before = sequence.substr(boundaries[0] - 1, start_loc - boundaries[0]);
        after = sequence.substr(end_loc, boundaries[1] - end_loc);
        variated_seq = before + after;
        break;
    default:
        break;
    }

    delete[] boundaries;

    return variated_seq;
}

/**
 * @brief Integrates a variation into the sequence graph by creating bubbles for differences.
 * @param original_seq The sequence graph representing the original sequence.
 * @param variated_seq The sequence after applying the variation.
 * @param start_loc The starting position for the variation in the sequence graph.
 * @param end_loc The ending position for the variation in the sequence graph.
 */
void variate(sequence_graph& original_seq, std::string variated_seq, int start_loc, int end_loc, std::vector<int> SN_ids) {

    int* boundaries = find_boundaries(start_loc, end_loc);
    lcp::lps* var_str = new lcp::lps(variated_seq, false);
    var_str->deepen(level);
  

    variated_strs.push_back(var_str);

    // set the curr as the starting core
    core_node* curr = original_seq.head;
    for (int i = 0; i < boundaries[2] - 1; i++) {
        curr = curr->next.at(0);
    }

    int first_difference, last_difference_org, last_difference_var;
    
    core_node* starting_core = align_variation(curr, var_str, boundaries[2], boundaries[3], 
        first_difference, last_difference_org, last_difference_var);
    curr = starting_core;

    std::cout << "fd: " << first_difference << " ld: " << last_difference_var << " ld2: " << last_difference_org << std::endl;


    // add all the different nodes as a bubble
    for (int i = (first_difference < 0) ? 0 : first_difference; i <= last_difference_var; i++) {
        core_node* new_cn = new core_node();
        new_cn->core_value = &var_str->cores->at(i);
        new_cn->SN_ids = SN_ids;
        if (first_difference < 0 && boundaries[2] == 0) {
            new_cn->id = "bh" + std::to_string(global_no_of_bubbles) + "-" + std::to_string(i - first_difference);
        } else {
            new_cn->id = "b" + std::to_string(global_no_of_bubbles) + "-" + std::to_string(i - first_difference);
        }
        new_cn->SO = new_cn->core_value->start;
        new_cn->rank = false;
        curr->next.push_back(new_cn);
        curr = new_cn;
    }

    global_no_of_bubbles++;

    // connect the end of the bubble to the node they are the same again
    core_node* end_of_bubble = curr;
    curr = starting_core->next.at(0);

    bool last_core = false;
    for (int i = first_difference; i < last_difference_org; i++) {
        if (curr->next.size() > 0)
            curr = curr->next.at(0);
        else
            last_core = true;

    }   
    if (!last_core) 
        end_of_bubble->next.push_back(curr);

    delete[] boundaries;
}

    void print_usage() {
        std::cout << "Usage: ./vg_lcp -f <fasta_file> -v <vcf_file> -r <rgfa_path> -l <LCP_level>\n\n";
    }

/**
 * @brief Main function to process the input, construct the sequence graph, and apply variations.
 * @param argc The number of command-line arguments.
 * @param argv The array of command-line arguments.
 * @return Exit status of the program (0 for success, non-zero for errors).
 */
int main(int argc, char* argv[]) {

    int opt;
    bool fasta = false, vcf = false, rgfa = false, level = false;

    while ((opt = getopt(argc, argv, "f:v:l:r:")) != -1) {
        switch (opt) {
        case 'f':
            fasta_path = optarg; // Handle file argument
            break;
        case 'v':
            vcf_path = optarg; // Handle VCF argument
            break;
        case 'r':
            rgfa_path = optarg;
            break;
        case 'l':
            level = std::stoi(optarg); // Convert level argument to integer
            break;
        case '?': // Unknown option or missing argument
            std::cerr << "Unknown or missing argument for option: " << static_cast<char>(optopt) << "\n";
            return 1;
        default:
            std::cerr << "Error parsing options\n";
            print_usage();
            return 1;
        }
    }

    if (!fasta || !vcf || !rgfa || !level) {
        print_usage();
        return 1;
    } 

    lcp::encoding::init();
    fasta_content fc;
    std::vector<variation*> variation_list;
    if (!read_fasta(fasta_path, fc)) return 1;
    if (!read_vcf(vcf_path, variation_list, chrmsms)) return 1;

    lcp::lps* str = new lcp::lps (fc.sequence);
    str->deepen(level);

    sequence_graph sg;

    core_node* prev = nullptr;
    std::vector<int> dum;
    dum.push_back(-1);

    int indx_incr = idx_size;

    for (size_t i = 0; i < str->size(); i++) {
        core_node* new_cn = new core_node();  
        new_cn->core_value = &str->cores->at(i);
        new_cn->id = "o" + std::to_string(i);
        new_cn->rank = true;
        new_cn->SN_ids = dum;
        new_cn->SO = new_cn->core_value->start;
        if (i == str->size() - 1) new_cn->end_flag = true;

        if (sg.head == nullptr) {
            sg.head = new_cn;
        }
        if (prev != nullptr) {
            prev->next.push_back(new_cn);
        }

	if ((int) new_cn->core_value->start >= curr_indx) {
	   indices.push_back(new_cn);
	   curr_indx += indx_incr;
	}

        prev = new_cn;
        sg.size++;
    }
    float count = 0;
    float  total_size = (float)(variation_list.size());	
    for (variation* v : variation_list) {
	std::cout << ((count / total_size) * 100) << " : " << count << " / " << total_size <<  std::endl;
	count++;
        if (strcmp(v->alt.c_str(), ".") == 0) { // deletion
      
            std::string variated_seq = 
                construct_variated_seq(fc.sequence, v->alt, v->pos, v->ref.size() + v->pos - 1, 2);
            variate(sg, variated_seq, v->pos, v->ref.size() + v->pos - 1, v->chromosom_ids);
        } else {
     
            if (strcmp(v->ref.c_str(), ".") == 0) { // insertion
                std::string variated_seq = 
                    construct_variated_seq(fc.sequence, v->alt, v->pos, v->pos, 0);
                variate(sg, variated_seq, v->pos, v->pos, v->chromosom_ids);
            } 
                
            else {
                std::string variated_seq = 
                    construct_variated_seq(fc.sequence, v->alt, v->pos, v->ref.size() + v->pos - 1, 1);
                variate(sg, variated_seq, v->pos, v->ref.size() + v->pos - 1, v->chromosom_ids);
            }
        }
    }

    print_seq(sg, rgfa_path, fc.chromosom);

    for (segment* seg : segments) {
        delete seg;
    }

    for (g_link* glnk : links) {
        delete glnk;
    }

    std::unordered_set<core_node*> visited;
    delete_core_nodes(sg.head, visited);

    for (variation* v : variation_list) {
        delete v;
    }

    for (lcp::lps* vs : variated_strs) {
        delete vs;
    }

    delete str;

    return 0;
}



