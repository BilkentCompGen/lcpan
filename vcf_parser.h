#ifndef VCF_PARSER
#define VCF_PARSER

#include <string>
#include <vector>

struct variation {
    int pos;
    std::string id;
    std::string ref;
    std::string alt;
    std::string chromosom;
    std::vector<int> chromosom_ids;
};

bool read_vcf(const std::string& file_name, std::vector<variation*>& variation_list, std::vector<std::string>& chrmsms);

#endif