#include <fstream>
#include <iostream>
#include <regex>
#include "bayes.hpp"


// Setup processing: load input files and define MPI task workload
void Bayes::setup_processing() {
    

}

// Check for minimal setup: a bed file + a dim file + phen file(s)
void Bayes::check_options() {
    
    std::cout << "will check passed options for completeness" << std::endl;
    
    if (opt.get_bed_file() == "") {
        std::cout << "FATAL  : no bed file provided! Please use the --bedfile option." << std::endl;
        exit(EXIT_FAILURE);
    }
    std::cout << "  bed file: OK - " << opt.get_bed_file() << "\n";

    if (opt.get_dim_file() == "") {
        std::cout << "FATAL  : no dim file provided! Please use the --dimfile option." << std::endl;
        exit(EXIT_FAILURE);
    }
    std::cout << "  dim file: OK - " << opt.get_dim_file() << "\n";

    if (opt.count_phen_files() == 0) {
        std::cout << "FATAL  : no phen file(s) provided! Please use the --phenfile option." << std::endl;
        exit(EXIT_FAILURE);
    }
    std::cout << "  phen file(s): OK - " << opt.count_phen_files() << " files passed.\n";
    opt.list_phen_files();
}


void Bayes::read_dim_file(int Nt, int Mt) {

    std::ifstream infile(opt.get_dim_file());
    std::string line;
    std::regex re("\\s+");

    if (infile.is_open()) {
        getline(infile, line);
        sregex_token_iterator first{line.begin(), line.end(), re, -1}, last;
        std::vector<std::string> tokens{first, last};
        infile.close();
        if (tokens.size() != 2) {
            std::cout << "FATAL: dim file should contain a single line with 2 integers" << std::endl;
            exit(EXIT_FAILURE);
        }
        Nt = atoi(tokens[0].c_str());
        Mt = atoi(tokens[1].c_str());
    } else {
        std::cout << "FATAL: could not open dim file: " << opt.get_dim_file() << std::endl;
        exit(EXIT_FAILURE);
    }
}
