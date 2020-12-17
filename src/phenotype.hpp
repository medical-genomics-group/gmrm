#pragma once
#include <string>
#include <vector>
#include <memory>
#include "options.hpp"
#include "distributions.hpp"
//#include "utilities.hpp"


class Phenotype {
    
public:
    Phenotype(std::string fp, const Options& opt, const int N, const int M);

    Phenotype(const Phenotype& rhs);

    ~Phenotype() {
        //std::cout << "\\\\\\ calling Phenotype dtor" << std::endl;
        if (mave    != nullptr)  _mm_free(mave);
        if (msig    != nullptr)  _mm_free(msig);
        if (epsilon != nullptr)  _mm_free(epsilon);
    }
    std::string get_filepath() const { return filepath; }
    void print_info() const;
    std::vector<unsigned char>& get_mask4() { return mask4; }
    std::vector<int>&           get_midx()  { return midx;  }
    std::vector<double>&        get_denom() { return denom; }
    std::vector<double>&        get_muk()   { return muk;   }
    std::vector<double>&        get_logl()  { return logl;  }
    std::vector<double>&        get_acum()  { return acum;  }
    std::vector<int>&           get_comp()  { return comp;  }
    std::vector<std::vector<double>>&  get_pi_est() { return pi_est; }
    std::vector<std::vector<int>>&     get_cass()   { return cass; }
    
    int     get_nas()   const { return nas; }
    int     get_nonas() const { return nonas; }
    double* get_mave()        { return mave; }
    double* get_msig()        { return msig; }
    double* get_epsilon()     { return epsilon; }
    double  get_epssum()      { return epssum; }
    double  get_sigmae()      { return sigmae; }
    double  get_sigmag()      { return sigmag; }
    double  get_mu()          { return mu; }

    void offset_epsilon(const double);
    void update_epsilon_sum();
    void update_epsilon_sigma();

    void set_rng(const unsigned int);
    void set_midx();
    void shuffle_midx();
    double sample_norm_rng();
    double sample_norm_rng(const double a, const double b);
    double sample_beta_rng();
    double sample_beta_rng(const double a, const double b);
    double sample_unif_rng();
    unsigned int get_random_int() { return dist.get_random_number(); }
        
    void set_pi_est(const std::vector<std::vector<double>> val) { pi_est = val; }
    void set_sigmag(const double val) { sigmag = val; }
    void set_mu(const double val) { mu = val; }

    void   set_marker_acum(const int idx, const double val) { acum[idx] = val; } 
    double get_marker_acum(const int idx) { return acum[idx]; } 

    void   set_marker_beta(const int idx, const double val) { betas[idx] = val; }
    double get_marker_beta(const int idx) { return betas[idx]; }

    int    get_marker_local_index(const int shuff_idx);
    double get_marker_ave(const int idx) { return mave[idx]; }
    double get_marker_sig(const int idx) { return msig[idx]; }

    void update_epsilon(const double* dbeta, const unsigned char* bed);

private:
    Distributions dist;
    std::string filepath;
    int nonas = 0;
    int nas   = 0;
    int im4   = 0;
    const unsigned int M = 0;
    const unsigned int N = 0;
    std::vector<double> betas;
    std::vector<double> data;
    std::vector<unsigned char> mask4;
    std::vector<int> midx;
    std::vector<double> denom;
    std::vector<double> muk;
    std::vector<double> logl;
    std::vector<double> acum;
    std::vector<int> comp;
    std::vector<std::vector<double>> pi_est;
    std::vector<std::vector<int>> cass;
    double* mave    = nullptr;
    double* msig    = nullptr;
    double* epsilon = nullptr; // starts with centered normalized phenotype
    double epssum = 0.0;
    double sigmae = 0.0;
    double sigmag = 0.0;
    double mu     = 0.0;
    void read_file(const Options& opt);
    friend class PhenMgr;
};


class PhenMgr {

public:
    PhenMgr() = default;
    PhenMgr(const Options& opt, const int N, const int M) {
        //std::cout << "+++Calling PhenMgr ctor" << std::endl; 
        read_phen_files(opt, N, M);
        //std::cout << "---PhenMgr ctor done" << std::endl; 
    }
    //~PhenMgr() {
    //    std::cout << "/!\\ Calling PhenMgr dtor" << std::endl; 
    //}
    void print_info();
    void compute_markers_statistics(const unsigned char* bed_data, const int N, const int M, const int mbytes);
    std::vector<Phenotype>& get_phens() { return phens; }
    void display_markers_statistics(const int n);
    void read_phen_files(const Options& opt, const int N, const int M);

private:
    std::vector<Phenotype> phens;
};

