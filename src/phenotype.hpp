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
    void epsilon_stats();
    void set_rng(const unsigned int);
    void set_midx();
    void shuffle_midx();
    void sample_mu_norm_rng();
    void sample_sigmag_beta_rng();
    void sample_sigmag_beta_rng(const double a, const double b);
    unsigned int get_random_int() { return dist.get_random_number(); }
    double get_beta(const int idx) { return betas[idx]; }
    int    get_marker_local_index(const int shuff_idx);
    double get_marker_ave(const int idx) { return mave[idx]; }
    double get_marker_sig(const int idx) { return msig[idx]; }

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

