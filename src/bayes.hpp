#pragma once
#include <iostream>
#include <typeinfo>
#include "options.hpp"
#include "phenotype.hpp"
#include "dimensions.hpp"


class Bayes {

public:

    const double V0E  = 0.0001;
    const double S02E = 0.0001;
    const double V0G  = 0.0001;
    const double S02G = 0.0001;


    Bayes(const Options& opt, const Dimensions& dims) : opt(opt),
                                                        rank(dims.get_rank()),
                                                        nranks(dims.get_nranks()),
                                                        N(dims.get_nt()),
                                                        Mt(dims.get_mt()),
                                                        K(opt.get_s().size() + 1),
                                                        G(opt.get_ngroups())  {

        cva.resize(G);
        cvai.resize(G);
        pi_prior.resize(G);
        mtotgrp.resize(G);

        set_block_of_markers();

        for (int i=0 ; i<G; i++) {
            mtotgrp.at(i) = 0;
            cva[i].resize(K, 0);
            cvai[i].resize(K, 0);
            double sum_cva = 0.0;
            for (int j=0; j<K-1; j++) {
                cva[i][j+1]  = opt.get_s().at(j);
                cvai[i][j+1] = 1.0 / cva[i][j+1];
                sum_cva += cva[i][j+1];
            }

            pi_prior[i].resize(K, 0.5);
            for (int j=1; j<K; j++) {
                pi_prior[i][j] =  pi_prior[i][0] * cva[i][j] / sum_cva;
                //printf("pi_prior[%d][%d] = %20.15f\n", i, j+1, pi_prior[i][j+1]);
            }
        }        

        setup_processing();
    }

    ~Bayes() {
        //std::cout << "## calling Bayes dtor" << std::endl;
        if (bed_data != nullptr)  _mm_free(bed_data);
    }

    void process();
    double dot_product(const int mloc, double* phen, const double mu, const double sigma);
    void list_phen_files() const { opt.list_phen_files(); }
    int  get_N()  { return N;  } // Invariant over tasks
    int  get_M()  { return M; } // Number of markers processed by task
    int  get_Mt() { return Mt; } // Total number of markers, sum over tasks
    int  get_Mm() { return Mm; } // Maximum number of markers per task (others may have M + 1)
    int  get_K()  { return K;  }
    void shuffle_markers();
    int  get_marker_group(const int mglob) { return 0; } //todo: adpat when groups are activated
    void update_epsilon(const int* counts, const double* dbetas, const unsigned char* recv_bed);
    void check_openmp();

private:
    const Options opt;
    PhenMgr pmgr;
    const int N = 0;
    const int Mt = 0;
    const int rank = 0;
    const int nranks = 0;
    unsigned char* bed_data = nullptr;
    const int K = 0;
    const int G = 0;
    
    std::vector<int> mtotgrp;
    std::vector<int> groups;
    std::vector<std::vector<double>> cva;
    std::vector<std::vector<double>> cvai;
    std::vector<std::vector<double>> pi_prior;

    int S = 0;              // task marker start 
    int M = 0;              // task marker length
    int Mm = 0;
    size_t mbytes = 0;

    void check_options();
    void setup_processing();
    void set_block_of_markers();
    void load_genotype();
    void check_processing_setup();
    void read_group_file();
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
