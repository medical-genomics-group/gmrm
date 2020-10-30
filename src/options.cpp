#include "options.hpp"
#include <iostream>
#include <cstring>
#include <sstream>
#include <iostream>
#include <fstream>

using namespace std;

// Function to parse command line options
void Options::read_command_line_options(int argc, char** argv) {

    std::stringstream ss;

    for (int i=1; i<argc; ++i) {

        // List of phenotype files to read; comma separated if more than one.
        if (!strcmp(argv[i], "--pheno")) {
            std::string phen_files = argv[++i];
            ss << "--pheno " << phen_files << "\n";
            std::stringstream ss(phen_files);
            std::string filepath;
            while (getline(ss, filepath, ',')) {
                ifstream phenfile(filepath);
                if (phenfile.is_open()) {
                    
                } else {
                    std::cout << "FATAL: file " << filepath << " not found\n";
                    exit(EXIT_FAILURE);
                }
                phenotypes.push_back(filepath);
            }
        }
    }

    std::cout << ss.str() << std::endl;
}

void Options::print_phenotypes(void) {
    for (auto f = phenotypes.begin(); f != phenotypes.end(); ++f) {
        std::cout << "input phenotype file: " << *f << std::endl;
    } 
}

