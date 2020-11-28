#pragma once

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/beta_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <boost/random/normal_distribution.hpp>

class Distributions {

public:

    void set_rng(unsigned int seed) {
        rng = boost::mt19937(seed);
    }

    unsigned int get_random_number() {
        return rng();
    }
 
    boost::mt19937& get_rng() {
        return rng;
    }

    double beta_rng(const double a, const double b) {
        //std::cout << "@@ beta_rng " << a << ", " << b << std::endl;
        boost::random::beta_distribution<double> mybeta(a, b);
        boost::random::variate_generator<boost::mt19937&, boost::random::beta_distribution<> > rand_beta(rng, mybeta);
        double val = rand_beta();
        //std::cout << "@@ beta_rng val = " << val << std::endl;
        return val;
    }

    double norm_rng(double mean, double sigma2) {
        //std::cout << "@@ norm_rng on " << mean << ", " << sigma2 << std::endl;
        boost::normal_distribution<double> nd(mean, std::sqrt(sigma2));
        boost::variate_generator< boost::mt19937&, boost::normal_distribution<> > var_nor(rng, nd);
        return var_nor();
    }

private:
    boost::mt19937 rng;
};


