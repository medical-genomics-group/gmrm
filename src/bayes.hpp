#pragma once
#include <iostream>
#include <typeinfo>
#include "options.hpp"

class Bayes {

public:
    Bayes() = default;
    Bayes(Options inopt) : opt(inopt) {
        check_options();
    }
    void list_phen_files() const { opt.list_phen_files(); }
    void setup_processing();

private:
    const Options opt;
    void check_options();
    void read_dim_file(int Nt, int Mt);
};


class BayesRR : public Bayes {

public:
    BayesRR(Options inopt) : Bayes(inopt) {
    }    
};


class BayesFH : public Bayes {

public:
    BayesFH(Options inopt) : Bayes(inopt) {
    }

};
