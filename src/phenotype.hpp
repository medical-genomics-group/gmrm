#pragma once
#include <string>
#include <vector>
#include <memory>
#include "options.hpp"


class Phenotype {
    
public:
    Phenotype(std::string fp, const Options& opt) : filepath(fp), nonas(0), nas(0) {
        //std::cout << ">>>>> calling Phenotype ctor on " << filepath << std::endl;
        read_file(opt);
        //std::cout << "<<<<< Phenotype ctor done." << filepath << std::endl;
    }

    Phenotype(const Phenotype& rhs);

    ~Phenotype() {
        //std::cout << "\\\\\\ calling Phenotype dtor" << std::endl;
        if (mave != nullptr)  _mm_free(mave);
        if (msig != nullptr)  _mm_free(msig);
        //std::cout << "????????????? csp = " << csp << " " << filepath << std::endl;
        if (csp  != nullptr)  _mm_free(csp);
        //std::cout << "///// calling Phenotype dtor" << std::endl;
    }
    std::string get_filepath() const { return filepath; }
    void print_info() const;
    std::vector<unsigned char>& get_mask4() { return mask4; };
    int get_nas()   const { return nas; }
    int get_nonas() const { return nonas; }
    double* get_mave() { return mave; }
    double* get_msig() { return msig; }
    double* get_csp()  { return csp; }
    void offset_epsilon(const double);

private:
    std::string filepath;
    int nonas = 0;
    int nas = 0;
    int im4 = 0;
    const unsigned int M = 0;
    std::vector<double> data;
    std::vector<unsigned char> mask4;
    double* mave = nullptr;
    double* msig = nullptr;
    double* csp  = nullptr; // csp: Centered and Scaled Phenotype
    std::vector<int> midx;

    void read_file(const Options& opt);

    friend class PhenMgr;
};


class PhenMgr {

public:
    PhenMgr() = default;
    PhenMgr(const Options& opt) {
        //std::cout << "+++Calling PhenMgr ctor" << std::endl; 
        read_phen_files(opt);
        //std::cout << "---PhenMgr ctor done" << std::endl; 
    }
    //~PhenMgr() {
    //    std::cout << "/!\\ Calling PhenMgr dtor" << std::endl; 
    //}
    void print_info();
    void compute_markers_statistics(const unsigned char* bed_data, const int N, const int M, const int mrk_bytes);
    std::vector<Phenotype> get_phens() { return phens; }
    void display_markers_statistics(const int n);
    void read_phen_files(const Options& opt);

private:
    std::vector<Phenotype> phens;
};

