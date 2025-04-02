#ifndef __UTILS_H__
#define __UTILS_H__

#include "struct_def.h"
#include "lps.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <pthread.h>

#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

#define SV_LEN_BOUNDARY 50

/**
 * Performs binary search on a sorted array to find the index of a given key.
 *
 * @param arr  Pointer to the sorted array.
 * @param size Number of elements in the array.
 * @param key  The value to search for.
 * @return Index of the key if found, otherwise -1.
 */
int binary_search(uint64_t *arr, uint64_t size, uint64_t key);

/**
 * Sorts an array using the quicksort algorithm.
 *
 * @param array Pointer to the array to be sorted.
 * @param low   Starting index of the array (or subarray).
 * @param high  Ending index of the array (or subarray).
 */
void quicksort(uint64_t *arr, int low, int high);

/**
 * Prints three sequences as single segment in GFA or rGFA format.
 *
 * @param id       Sequence identifier.
 * @param seq1     The first nucleotide sequence.
 * @param seq1_len Length of the first sequence.
 * @param seq2     The second nucleotide sequence.
 * @param seq2_len Length of the second sequence.
 * @param seq3     The third nucleotide sequence.
 * @param seq3_len Length of the third sequence.
 * @param seq_name Name of the sequence.
 * @param start    Start position of the sequence.
 * @param rank     Rank of the sequence.
 * @param is_rgfa  Flag to determine if rGFA format should be used.
 * @param out      Output file stream.
 */
void print_seq3(uint64_t id, const char *seq1, int seq1_len, const char *seq2, int seq2_len, const char *seq3, int seq3_len, const char *seq_name, int start, int rank, int is_rgfa, FILE *out);

/**
 * Prints two sequences as single segment in GFA or rGFA format.
 *
 * @param id       Sequence identifier.
 * @param seq1     The first nucleotide sequence.
 * @param seq1_len Length of the first sequence.
 * @param seq2     The second nucleotide sequence.
 * @param seq2_len Length of the second sequence.
 * @param seq_name Name of the sequence.
 * @param start    Start position of the sequence.
 * @param rank     Rank of the sequence.
 * @param is_rgfa  Flag to determine if rGFA format should be used.
 * @param out      Output file stream.
 */
void print_seq2(uint64_t id, const char *seq1, int seq1_len, const char *seq2, int seq2_len, const char *seq_name, int start, int rank, int is_rgfa, FILE *out);

/**
 * Prints a sequence in GFA or rGFA format.
 *
 * @param id       Sequence identifier.
 * @param seq      The nucleotide sequence.
 * @param seq_len  Length of the sequence.
 * @param seq_name Name of the sequence.
 * @param start    Start position of the sequence.
 * @param rank     Rank of the sequence.
 * @param is_rgfa  Flag to determine if rGFA format should be used.
 * @param out      Output file stream.
 */
void print_seq(uint64_t id, const char *seq, int seq_len, const char *seq_name, int start, int rank, int is_rgfa, FILE *out);

/**
 * Prints three sequences as single segment in GFA or rGFA format. Unlike `print_seq`, this function
 * does not print index, but prints >id to index information (`SN:Z:`)
 *
 * @param id       Sequence identifier.
 * @param seq1     The first nucleotide sequence.
 * @param seq1_len Length of the first sequence.
 * @param seq2     The second nucleotide sequence.
 * @param seq2_len Length of the second sequence.
 * @param seq3     The third nucleotide sequence.
 * @param seq3_len Length of the third sequence.
 * @param seq_name The name/ID of the sequence.
 * @param order    The index of the sequence (order).
 * @param start    Start position of the sequence.
 * @param rank     Rank of the sequence.
 * @param is_rgfa  Flag to determine if rGFA format should be used.
 * @param out      Output file stream.
 */
void print_seq3_vg(uint64_t id, const char *seq1, int seq1_len, const char *seq2, int seq2_len, const char *seq3, int seq3_len, const char *seq_name, int order, int start, int rank, int is_rgfa, FILE *out);

/**
 * Prints two sequences as single segment in GFA or rGFA format. Unlike `print_seq`, this function
 * does not print index, but prints >id to index information (`SN:Z:`)
 *
 * @param id       Sequence identifier.
 * @param seq1     The first nucleotide sequence.
 * @param seq1_len Length of the first sequence.
 * @param seq2     The second nucleotide sequence.
 * @param seq2_len Length of the second sequence.
 * @param seq_name The name/ID of the sequence.
 * @param order    The index of the sequence (order).
 * @param start    Start position of the sequence.
 * @param rank     Rank of the sequence.
 * @param is_rgfa  Flag to determine if rGFA format should be used.
 * @param out      Output file stream.
 */
void print_seq2_vg(uint64_t id, const char *seq1, int seq1_len, const char *seq2, int seq2_len, const char *seq_name, int order, int start, int rank, int is_rgfa, FILE *out);

/**
 * Prints a sequence in GFA or rGFA format. Unlike `print_seq`, this function
 * does not print index, but prints >id to index information (`SN:Z:`)
 *
 * @param id       Sequence identifier.
 * @param seq      The nucleotide sequence.
 * @param seq_len  Length of the sequence.
 * @param seq_name The name/ID of the sequence.
 * @param order    The index of the sequence (order).
 * @param start    Start position of the sequence.
 * @param rank     Rank of the sequence.
 * @param is_rgfa  Flag to determine if rGFA format should be used.
 * @param out      Output file stream.
 */
void print_seq_vg(uint64_t id, const char *seq, int seq_len, const char *seq_name, int order, int start, int rank, int is_rgfa, FILE *out);

/**
 * Prints a link between two sequences in GFA format.
 *
 * @param id1        Identifier of the first sequence.
 * @param sign1      Orientation of the first sequence ('+' or '-').
 * @param id2        Identifier of the second sequence.
 * @param sign2      Orientation of the second sequence ('+' or '-').
 * @param overlap    Length of the overlap between the sequences.
 * @param out        Output file stream.
 */
void print_link(uint64_t id1, char sign1, uint64_t id2, char sign2, uint64_t overlap, FILE *out);

/**
 * Finds the latest core index before a given range and the first core index after it.
 *
 * @param start_loc         Start location of the range.
 * @param end_loc           End location of the range.
 * @param chrom             Pointer to the chromosome structure containing core regions.
 * @param start_index       Estimated start index of LCP core to search.
 * @param latest_core_index Pointer to store the latest core index before start_loc.
 * @param first_core_after  Pointer to store the first core index after end_loc.
 */
void find_boundaries(uint64_t start_loc, uint64_t end_loc, const struct chr *chrom, uint64_t start_index, uint64_t *latest_core_index, uint64_t *first_core_after);

/**
 * Modifies start indices of the LCP cores if no overlap is allowed.
 * 
 * @param seqs       The LCP cores to be refined.
 * @param no_overlap Boolean indicator that makes refinement is it is set to 1.
 */
void refine_seqs(struct ref_seq *seqs, int no_overlap);

/**
 * Modifies start indices of the LCP cores if no overlap is allowed.
 * 
 * @param str        The LCP cores to be refined.
 * @param no_overlap Boolean indicator that makes refinement is it is set to 1.
 */
void refine_seq(struct lps *str, int no_overlap);

/**
 * Prints all Paths in given sequences. The ids should be initialized
 * 
 * @param ref_seq   The reference sequences
 * @param out       Output file to write path.
 */
void print_path(struct ref_seq *seqs, FILE *out);

#endif