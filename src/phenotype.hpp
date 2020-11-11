#pragma once
#include <string>
#include <vector>
#include <memory>
#include "options.hpp"

class Phenotype {
    
public:
    Phenotype(std::string fp, const Options& opt) : filepath(fp), nonas(0), nas(0) {
        read_file(opt);
    }
    ~Phenotype() {
        if (mavg != nullptr)  _mm_free(mavg);
        if (mstd != nullptr)  _mm_free(mstd);
    }
    std::string get_filepath() const { return filepath; }
    void print_info() const;
    std::vector<unsigned char>& get_mask8() { return mask8; }; 
    std::vector<unsigned char>& get_mask4() { return mask4; };
    int get_nas()   const { return nas; }
    int get_nonas() const { return nonas; }

private:
    std::string filepath;
    int nonas;
    int nas;
    std::vector<double> data;
    std::vector<unsigned char> mask8; // 8 inds per byte
    std::vector<unsigned char> mask4; // 4 inds per byte
    double* mavg = nullptr;
    double* mstd = nullptr;
    void read_file(const Options& opt);

    friend class PhenMgr;
};


class PhenMgr {

public:
    PhenMgr(const Options& opt) {
        read_phen_files(opt);
    }
    void print_info();
    void compute_markers_statistics(const unsigned char* bed_data, const int N, const int M, const int mrk_bytes);
    std::vector<Phenotype>& get_phens() { return phens; }

private:
    std::vector<Phenotype> phens;
    void read_phen_files(const Options& opt);
};

