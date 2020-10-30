#include <iostream>
#include <fstream>
#include <regex>
#include <cmath>
#include "phenotype.hpp"

using namespace std;


// Read phenotype file
// Assume PLINK format: Family ID, Individual ID, Phenotype
// One row per individual
void Phenotype::read_file() {

    ifstream infile(filepath);
    string line;
    regex re("\\s+");

    if (infile.is_open()) {
        nonas = 0, nas = 0;
        while (getline(infile, line)) {
            sregex_token_iterator first{line.begin(), line.end(), re, -1}, last;
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
        cout << "FATAL: could not open phenotype file: " << filepath << endl;
        exit(EXIT_FAILURE);
    }

    //assert(nonas + nas == numInds);
}

void Phenotype::print_info() const {
    printf("INFO   : %s has %d NAs and %d non-NAs.\n", get_filepath().c_str(), nas, nonas);
}
