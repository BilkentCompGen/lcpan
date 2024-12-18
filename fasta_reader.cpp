#include <iostream>
#include <fstream>

#include "fasta_reader.h"

bool read_fasta(const std::string& file_name, fasta_content& fasta_c) {
    std::ifstream fasta_file(file_name);

    if (!fasta_file.is_open()) {
        std::cerr << "Failed to open the file!" << std::endl;
        return false;
    }

    std::string line;
    
    while (std::getline(fasta_file, line)) {
        if (line[0] == '>') {
            fasta_c.chromosom = line.substr(1, line.size() - 1);
        } else {
            fasta_c.sequence += line;
        }
    }


    return true;
}
