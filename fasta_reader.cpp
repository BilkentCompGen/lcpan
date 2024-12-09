#include <iostream>
#include <fstream>
#include <string>

bool read_fasta(const std::string& file_name, std::string& sequence) {
    std::ifstream fasta_file(file_name);

    if (!fasta_file.is_open()) {
        std::cerr << "Failed to open the file!" << std::endl;
        return false;
    }

    std::string line;
    
    while (std::getline(fasta_file, line)) {
        if (line[0] == '>') {
            std::cout << "Header: " << line << std::endl;  
        } else {
            sequence += line;
        }
    }

    return true;
}
