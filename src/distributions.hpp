#pragma once

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/variate_generator.hpp>


class Distributions {

public:
    Distributions() {
    }
    void set_rng(unsigned int seed) { rng = boost::mt19937(seed); }
    boost::mt19937& get_rng() { return rng; }

private:
    boost::mt19937 rng;
};


