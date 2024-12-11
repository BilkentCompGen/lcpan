#include <iostream>
#include <fstream>
#include <sstream>

#include "vcf_parser.h"

bool read_vcf(const std::string& file_name, std::vector<variation*>& variation_list) {
    std::ifstream vcf_file(file_name);

    if (!vcf_file.is_open()) {
        std::cerr << "Failed to open the file!";
        return false;
    }

    std::string line;
    while (std::getline(vcf_file, line)) {
        if (line.empty() || line[0] == '#') continue;

        variation* new_var = new variation;

        std::stringstream ss(line);
        std::string word;
        for (int i = 0; i < 5; ++i) {
            ss >> word;
            switch (i) {
            case 0:
                new_var->chromosom = word;
                break;
            case 1:
                new_var->pos = std::stoi(word);
                break;
            case 2:
                new_var->id = word;
                break;
            case 3:
                new_var->ref = word;
                break;
            case 4:
                new_var->alt = word;
                break;
            default:
                break;
            }
        }

        variation_list.push_back(new_var);
    }

    return true;
}