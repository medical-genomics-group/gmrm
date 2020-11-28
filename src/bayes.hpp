#pragma once
#include <iostream>
#include <typeinfo>
#include "options.hpp"
#include "phenotype.hpp"
#include "dimensions.hpp"


class Bayes {

public:
    Bayes(const Options&    opt,
          const Dimensions& dims) : opt(opt),
                                    rank(dims.get_rank()),
                                    nranks(dims.get_nranks()),
                                    N(dims.get_nt()),
                                    Mt(dims.get_mt()) {
        set_block_of_markers();
        setup_processing();
    }

    ~Bayes() {
        //std::cout << "## calling Bayes dtor" << std::endl;
        if (bed_data != nullptr)  _mm_free(bed_data);
    }

    void process();
    void list_phen_files() const { opt.list_phen_files(); }
    int  get_N()  { return N;  } // Invariant over tasks
    int  get_M()  { return M;  } // Number of markers processed by task
    int  get_Mt() { return Mt; } // Total number of markers, sum over tasks
    void shuffle_markers();


private:
    const Options opt;
    PhenMgr pmgr;
    const int N = 0;
    const int Mt = 0;
    const int rank = 0;
    const int nranks = 0;
    unsigned char* bed_data = nullptr;

    int S = 0;              // task marker start 
    int M = 0;              // task marker length
    size_t mrk_bytes = 0;
    size_t mrk_uints = 0;

    void check_options();
    void setup_processing();
    void set_block_of_markers();
    void load_genotype();
    void check_processing_setup();
};


class BayesRR : public Bayes {

public:
    
    BayesRR(const Options& opt, 
            const Dimensions& dims) : Bayes(opt,
                                            dims) { 
        //std::cout << "calling BayesRR constructor" << std::endl;
    }
};

/*
class BayesFH : public Bayes {

public:
BayesFH(const Options& opt, const PhenMgr& pmgr, const Dimensions& dims) : Bayes(opt, pmgr, rank, nranks) {
    }

};
*/
