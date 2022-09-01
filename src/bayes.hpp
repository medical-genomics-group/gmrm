#pragma once
#include <iostream>
#include <typeinfo>
#include <map>
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
                                                        //K(opt.get_nmixtures()),
                                                        G(opt.get_ngroups()) 
                                                        //cva(opt.get_cva()), 
                                                        //cvai(opt.get_cvai()) 
                                                        {

        std::cout << "Calling Bayes constructor" << std::endl;
        set_block_of_markers();

        if (!opt.predict()) {
            mtotgrp.resize(G);

            for (int i = 0; i < G; i++) {
                mtotgrp.at(i) = 0;
                //double sum_cva = 0.0;
                //sum_cva += cva[i];

            };
            pi_prior.resize(M);

            for (int i = 0; i < M; i++) {
                // Multitrate
                pi_prior[i] = 0.5;
            };
        }

        setup_processing();
        cout << "Setup processing done." <<"\n";
    }

    ~Bayes() {
        if (bed_data != nullptr)  _mm_free(bed_data);
    }

    //void predict();
    void process();
    void cross_bim_files();
    double dot_product(const int mloc, double* eps, const double mu, const double sigma);
    std::vector<std::vector<double>> eps_eps_product(const int q);
    std::vector<std::vector<double>> eps_eps_product1(const int q);
    void list_phen_files() const { opt.list_phen_files(); }
    int  get_N()  { return N;  } // Invariant over tasks
    int  get_M()  { return M;  } // Number of markers processed by task
    int  get_Mt() { return Mt; } // Total number of markers, sum over tasks
    int  get_Mm() { return Mm; } // Maximum number of markers per task (others may have M + 1)
    int  get_K()  { return K;  }
    //void shuffle_markers();

    int  get_marker_group(const int mglob) { 
        return group_index.at(mglob);
    }

    //void update_epsilon(const int* counts, const double* dbetas, const unsigned char* recv_bed);
    void update_epsilon(const int mloc, const std::vector<double> dbetas);
    std::vector<std::vector<double>> create_identity(const int q);
    void check_openmp();
    //void print_cva();
    //void print_cvai();

private:
    const Options opt;
    PhenMgr pmgr;
    const int N = 0;
    const int Mt = 0;
    const int rank = 0;
    const int nranks = 1;
    unsigned char* bed_data = nullptr;
    const int K = 0;
    const int G = 0;
    
    std::vector<int> mtotgrp;
    std::vector<int> group_index;

    // Multitrait
    std::vector<double> pi_prior; // Vector of priors for markers
    //const std::vector<double> cva; // Vector for groups
    //const std::vector<double> cvai; // Vector for groups

    std::vector<std::string>   rsid;
    std::map<std::string, int> m_refrsid;

    int S = 0;              // task marker start 
    int M = 0;              // task marker length
    int Mm = 0;
    size_t mbytes = 0;

    //void check_options();
    void setup_processing();
    void set_block_of_markers();
    void load_genotype();
    void check_processing_setup();
    void read_group_index_file(const std::string& file);
};


class BayesRR : public Bayes {

public:
    
    BayesRR(const Options& opt,const Dimensions& dims) : Bayes(opt,dims) { 
        std::cout << "Calling BayesRR constructor" << std::endl;
    }
};
