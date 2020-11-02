#pragma once
#include <iostream>
#include <typeinfo>
#include "options.hpp"
#include "phenotype.hpp"
#include "dimensions.hpp"


class Bayes {

public:
    //Bayes() = default;
    Bayes(const Options& opt, const PhenMgr& pmgr, const Dimensions& dims) : opt(opt), rank(dims.get_rank()), nranks(dims.get_nranks()), Nt(dims.get_nt()), Mt(dims.get_mt()) {
        std::cout << "calling Bayes constructor" << std::endl;
        check_options();
        setup_processing();
    }

    ~Bayes() {
        //std::cout << "rank " << rank << " calling Bayes destructor" << std::endl;
        if (bed_data != nullptr)  _mm_free(bed_data);
    }

    void list_phen_files() const { opt.list_phen_files(); }

private:
    const Options opt;
    const int Nt = 0;
    const int Mt = 0;
    const int rank = 0;
    const int nranks = 0;
    char* bed_data = nullptr;

    int S = 0;              // task marker start 
    int M = 0;              // task marker length
    size_t mrk_bytes = 0;
    size_t mrk_uints = 0;

    void check_options();
    void setup_processing();
    void set_block_of_markers();
    void load_genotype();
};


class BayesRR : public Bayes {

public:

    BayesRR(const Options& opt, const PhenMgr& pmgr, const Dimensions& dims) : Bayes(opt, pmgr, dims) { 
         std::cout << "calling BayesRR constructor" << std::endl;
    }

    /*
    BayesRR(const Options& opt, const PhenMgr& pmgr, const Dimensions& dims) : Bayes(opt, pmgr, rank, nranks) {
    } 
    */
};

/*
class BayesFH : public Bayes {

public:
BayesFH(const Options& opt, const PhenMgr& pmgr, const Dimensions& dims) : Bayes(opt, pmgr, rank, nranks) {
    }

};
*/
