#ifndef FASTA_READER
#define FASTA_READER

#include <string>

struct fasta_content {
    std::string sequence;
    std::string chromosom;
};



bool read_fasta(const std::string& file_name, fasta_content& fasta_c);

#endif