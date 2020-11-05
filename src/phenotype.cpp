#include <iostream>
#include <fstream>
#include <regex>
#include <cmath>
#include "phenotype.hpp"


void PhenMgr::print_info() const {
    for (auto phen = phen_mgr.begin(); phen != phen_mgr.end(); ++phen) {
        (*phen).print_info();
    }
}

void PhenMgr::read_phen_files(const Options& opt) {
    std::vector<std::string> phen_files = opt.get_phen_files();
    for (auto fp = phen_files.begin(); fp != phen_files.end(); ++fp) {
        if (opt.verbosity_level(3))
            std::cout << "Reading phenotype file: " << *fp << std::endl;
        Phenotype phenotype(*fp, opt);
        phen_mgr.push_back(phenotype);
    }
}


// Read phenotype file
// Assume PLINK format: Family ID, Individual ID, Phenotype
// One row per individual
void Phenotype::read_file(const Options& opt) {

    std::ifstream infile(filepath);
    std::string line;
    std::regex re("\\s+");

    if (infile.is_open()) {
        int line_n = 0;
        nonas = 0, nas = 0;
        while (getline(infile, line)) {            
            int m8 = line_n % 8;
            if (m8 == 0)  mask.push_back(0xFF);
            std::sregex_token_iterator first{line.begin(), line.end(), re, -1}, last;
            std::vector<std::string> tokens{first, last};
            if (tokens[2] == "NA") {
                nas += 1;
                data.push_back(NAN);
                if (opt.verbosity_level(3))
                    std::cout << " ... found NA on line " << line_n << ", m8 = " << m8 << " on byte " << int(line_n / 8) << std::endl;
                mask.at(int(line_n / 8)) &= ~(0b1 << m8);
            } else {
                nonas += 1;
                data.push_back(atof(tokens[2].c_str()));
            }

            line_n += 1;

            if (opt.verbosity_level(3)) {
                // Print fully handled previous byte
                if (m8 == 0 && line_n > 7 && line_n < 30) {
                    std::cout << "byte with NA " << unsigned(mask.at(int(line_n / 8) - 1)) << std::endl;
                }
            }

        }
        infile.close();
    } else {
        std::cout << "FATAL: could not open phenotype file: " << filepath << std::endl;
        exit(EXIT_FAILURE);
    }


    //assert(nonas + nas == numInds);
    
}

void Phenotype::print_info() const {
    printf("INFO   : %s has %d NAs and %d non-NAs.\n", get_filepath().c_str(), nas, nonas);
}
