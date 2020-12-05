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
                                    Mt(dims.get_mt()),
                                    K(opt.get_s().size() + 1),
                                    ngroups(opt.get_ngroups()) {
        set_block_of_markers();
        setup_processing();
        
        cva.resize(ngroups);
        cvai.resize(ngroups);
        for (int i=0 ; i<ngroups; i++) {
            cva[i].resize(K, 0);
            cvai[i].resize(K, 0);
            for (int j=0; j<K-1; j++) {
                cva[i][j+1]  = opt.get_s().at(j);
                cvai[i][j+1] = 1.0 / cva[i][j+1]; 
            }
        }
        std::cout << "cvai size = " << cvai.size() << std::endl;
        std::cout << "cvai size = " << cvai[0].size() << std::endl;
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
    int  get_Mm() { return Mm; } // Maximum number of markers per task (others may have M + 1)
    int  get_K()  { return K;  }
    void shuffle_markers();
    int  get_marker_group(const int mglob) { return 0; } //todo: adpat when groups are activated

private:
    const Options opt;
    PhenMgr pmgr;
    const int N = 0;
    const int Mt = 0;
    const int rank = 0;
    const int nranks = 0;
    unsigned char* bed_data = nullptr;
    const int K = 0;
    const int ngroups = 0;

    std::vector<std::vector<double>> cva;
    std::vector<std::vector<double>> cvai;

    int S  = 0;              // task marker start 
    int M  = 0;              // task marker length
    int Mm = 0;
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
