#pragma once
#include <string>
#include <vector>
#include "options.hpp"

class Phenotype {
    
private:
    std::string filepath;
    int nonas;
    int nas;
    std::vector<double> data;
    std::vector<unsigned char> mask;
    void read_file(const Options& opt);

public:
    Phenotype(std::string fp, const Options& opt) : filepath(fp), nonas(0), nas(0) {
        read_file(opt);
    }
    std::string get_filepath() const { return filepath; }
    void print_info() const;
};


class PhenMgr {

public:
    PhenMgr(const Options& opt) {
        read_phen_files(opt);
    }
    void print_info() const;

private:
    std::vector<Phenotype> phen_mgr;
    void read_phen_files(const Options& opt);
};

