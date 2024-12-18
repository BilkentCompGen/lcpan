#include <iostream>
#include <fstream>
#include <sstream>

#include "vcf_parser.h"

const char* ws = " \t\n\r\f\v";

// trim from end of string (right)
inline std::string& rtrim(std::string& s, const char* t = ws)
{
    s.erase(s.find_last_not_of(t) + 1);
    return s;
}

// trim from beginning of string (left)
inline std::string& ltrim(std::string& s, const char* t = ws)
{
    s.erase(0, s.find_first_not_of(t));
    return s;
}

// trim from both ends of string (right then left)
inline std::string& trim(std::string& s, const char* t = ws)
{
    return ltrim(rtrim(s, t), t);
}

bool read_vcf(const std::string& file_name, std::vector<variation*>& variation_list, std::vector<std::string>& chrmsms) {
    std::ifstream vcf_file(file_name);

    if (!vcf_file.is_open()) {
        std::cerr << "Failed to open the file!";
        return false;
    }

    std::string line;
    
    while (std::getline(vcf_file, line)) {

        if (line.empty() || (line[0] == '#' && line[1] == '#')) continue;

        if (line[0] == '#') {
            std::vector<std::string> tokens;

            std::istringstream iss(line);
            std::string token;

            while (std::getline(iss, token, ' ')) {
                tokens.push_back(trim(token));
            }

            bool format_seen = false;
            for (std::string t : tokens) {
                if (format_seen) {
                    chrmsms.push_back(t);
                }
                if (t.compare("FORMAT") == 0) {
                    format_seen = true;
                }
            }
            continue;
        }

        std::stringstream ss(line);
        std::string word;
        std::vector<variation*> v_list;
        variation* new_var = new variation;

        for (int i = 0; i < 8 + chrmsms.size(); ++i) {
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

            if (i > 8) {
                if (word[0] != 0 || word[3] != 0) {
                    new_var->chromosom_ids.push_back(i - 9);
                }
            }
        }

        std::vector<std::string> tokens;

        std::istringstream iss(new_var->alt);
        std::string token;

        while (std::getline(iss, token, ',')) {
            tokens.push_back(trim(token));
        }

        if (tokens.size() > 1) {
            for (std::string t : tokens) {
                variation* v = new variation; 
                v->alt = t;
                v->chromosom = new_var->chromosom;
                v->chromosom_ids = new_var->chromosom_ids;
                v->id = new_var->id;
                v->pos = new_var->pos;
                v->ref = new_var->ref;
                variation_list.push_back(v);
            }
        } else variation_list.push_back(new_var);
    }

    return true;
}