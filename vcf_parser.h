#ifndef VCF_PARSER
#define VCF_PARSER

#include <string>
#include <vector>

struct variation {
    std::string chromosom;
    int pos;
    std::string id;
    std::string ref;
    std::string alt;
};

bool read_vcf(const std::string& file_name, std::vector<variation*>& variation_list);

#endif