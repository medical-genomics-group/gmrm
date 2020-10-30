#pragma once

#include <string>
#include <vector>

using namespace std;


class Phenotype {
    
private:
    string filepath;
    int nonas;
    int nas;
    vector<double> data;
    void read_file();

public:
    Phenotype() = default;
    Phenotype(string fp) : filepath(fp), nonas(0), nas(0) {
        read_file();
    }
    string get_filepath() const { return filepath; }
    void print_info() const;
};
