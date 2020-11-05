#include <iostream>
#include <cstring>
#include <sstream>
#include <iostream>
#include <fstream>
#include "options.hpp"


// Function to parse command line options
void Options::read_command_line_options(int argc, char** argv) {

    std::stringstream ss;
    ss << "\nardyh command line options:\n";

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
        } else if (!strcmp(argv[i], "--verbosity")) {
            if (i == argc - 1) fail_if_last(argv, i);
            verbosity = atoi(argv[++i]);            
            ss << "--verbosity " << verbosity << "\n";
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

// Check for minimal setup: a bed file + a dim file + phen file(s)
void Options::check_options() {
    
    std::cout << "will check passed options for completeness" << std::endl;
    
    if (get_bed_file() == "") {
        std::cout << "FATAL  : no bed file provided! Please use the --bedfile option." << std::endl;
        exit(EXIT_FAILURE);
    }
    std::cout << "  bed file: OK - " << get_bed_file() << "\n";

    if (get_dim_file() == "") {
        std::cout << "FATAL  : no dim file provided! Please use the --dimfile option." << std::endl;
        exit(EXIT_FAILURE);
    }
    std::cout << "  dim file: OK - " << get_dim_file() << "\n";

    if (count_phen_files() == 0) {
        std::cout << "FATAL  : no phen file(s) provided! Please use the --phenfile option." << std::endl;
        exit(EXIT_FAILURE);
    }
    std::cout << "  phen file(s): OK - " << count_phen_files() << " files passed.\n";
    list_phen_files();
}
