#pragma once

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/beta_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>

class Distributions {

public:
    Distributions() {
    }
    void set_rng(unsigned int seed) { rng = boost::mt19937(seed); }
    boost::mt19937& get_rng() { return rng; }

    double beta_rng(double a, double b) {
        boost::random::beta_distribution<double> mybeta(a, b);
        boost::random::variate_generator<boost::mt19937&, boost::random::beta_distribution<> > rand_beta(rng, mybeta);
        return rand_beta();
    }

    double norm_rng(double mean, double sigma2) {
        boost::normal_distribution<double> nd(mean, std::sqrt(sigma2));
        boost::variate_generator<boost::mt19937&,boost::normal_distribution<> > var_nor(rng, nd);
        return var_nor();
    }

private:
    boost::mt19937 rng;
};


