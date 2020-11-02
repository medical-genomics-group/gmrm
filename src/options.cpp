#include <iostream>
#include <cstring>
#include <sstream>
#include <iostream>
#include <fstream>
#include "options.hpp"


// Function to parse command line options
void Options::read_command_line_options(int argc, char** argv) {

    std::stringstream ss;

    for (int i=1; i<argc; ++i) {

        if (!strcmp(argv[i], "--bedfile")) {
            if (i == argc - 1) fail_if_last(argv, i);
            bed_file = argv[++i];
            ss << "--bedfile " << bed_file << "\n";
        }
        else if (!strcmp(argv[i], "--dimfile")) {
            if (i == argc - 1) fail_if_last(argv, i);
            dim_file = argv[++i];
            ss << "--dimfile " << dim_file << "\n";
        }
        // List of phenotype files to read; comma separated if more than one.
        else if (!strcmp(argv[i], "--phenfiles")) {
            if (i == argc - 1) fail_if_last(argv, i);
            std::string cslist = argv[++i];
            ss << "--phenfiles " << cslist << "\n";
            std::stringstream sslist(cslist);
            std::string filepath;
            while (getline(sslist, filepath, ',')) {
                std::ifstream phen_file(filepath);
                if (phen_file.is_open()) {
                    phen_file.close();
                    phen_files.push_back(filepath);
                } else {
                    std::cout << "FATAL: file " << filepath << " not found\n";
                    exit(EXIT_FAILURE);
                }
            }
        } else {
            std::cout << "FATAL: option \"" << argv[i] << "\" unknown\n";
            exit(EXIT_FAILURE);
        }
    }

    std::cout << ss.str() << std::endl;
}

void Options::list_phen_files() const {
    for (auto phen = phen_files.begin(); phen != phen_files.end(); ++phen) {
        std::cout << " phen file: " << *phen << std::endl;
    } 
}

// Catch missing argument on last passed option
void Options::fail_if_last(char** argv, const int i) {
    std::cout << "FATAL  : missing argument for last option \"" << argv[i] <<"\". Please check your input and relaunch." << std::endl;
    exit(EXIT_FAILURE);
}
