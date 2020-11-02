#include <iostream>
#include <fstream>
#include <regex>
#include <cmath>
#include "phenotype.hpp"

void PhenMgr::read_phen_files(const Options& opt) {

    std::vector<std::string> phen_files = opt.get_phen_files();

    for (auto fp = phen_files.begin(); fp != phen_files.end(); ++fp) {
        std::cout << "will read phen " << *fp << std::endl;
        Phenotype phenotype(*fp);
        phen_mgr.push_back(phenotype);
    }
}


// Read phenotype file
// Assume PLINK format: Family ID, Individual ID, Phenotype
// One row per individual
void Phenotype::read_file() {

    std::ifstream infile(filepath);
    std::string line;
    std::regex re("\\s+");

    if (infile.is_open()) {
        nonas = 0, nas = 0;
        while (getline(infile, line)) {
            std::sregex_token_iterator first{line.begin(), line.end(), re, -1}, last;
            std::vector<std::string> tokens{first, last};
            if (tokens[2] == "NA") {
                nas += 1;
                data.push_back(NAN);
            } else {
                nonas += 1;
                data.push_back(atof(tokens[2].c_str()));
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
