#include <iostream>
#include <fstream>
#include <regex>
#include <cmath>
#include <immintrin.h>
#include <mpi.h>
#include "utilities.hpp"
#include "phenotype.hpp"
#include "dotp_lut.hpp"
#include "na_lut.hpp"
#include <boost/range/algorithm.hpp>
#include <boost/random/uniform_int.hpp>


Phenotype::Phenotype(std::string fp, const Options& opt, const int N, const int M) :
    filepath(fp),
    N(N),
    M(M),
    im4(N%4 == 0 ? N/4 : N/4+1),
    K(opt.get_nmixtures()),
    G(opt.get_ngroups()) {

    mave = (double*) _mm_malloc(size_t(M) * sizeof(double), 64);
    check_malloc(mave, __LINE__, __FILE__);
    msig = (double*) _mm_malloc(size_t(M) * sizeof(double), 64);
    check_malloc(msig, __LINE__, __FILE__);
    epsilon = (double*) _mm_malloc(size_t(im4*4) * sizeof(double), 64);
    check_malloc(epsilon, __LINE__, __FILE__);
    cass = (int*) _mm_malloc(G * K * sizeof(int), 64);
    check_malloc(cass, __LINE__, __FILE__);
    
    betas.assign(M, 0.0);
    acum.assign(M, 0.0);
    muk.assign(K, 0.0);
    denom.resize(K - 1);
    logl.resize(K);
    comp.assign(M, 0);
    beta_sqn.resize(G);
    m0.resize(G);
    sigmag.resize(G);

    dirich.clear();
    for (int i=0; i<K; i++) dirich.push_back(1.0);

    read_file(opt);
}

//copy ctor
Phenotype::Phenotype(const Phenotype& rhs) :
    dist(rhs.dist),
    filepath(rhs.filepath),
    nonas(rhs.nonas),
    nas(rhs.nas),
    im4(rhs.im4),
    N(rhs.N),
    M(rhs.M),
    G(rhs.G),
    K(rhs.K),
    betas(rhs.betas),
    data(rhs.data),
    mask4(rhs.mask4),
    midx(rhs.midx),
    denom(rhs.denom),
    muk(rhs.muk),
    logl(rhs.logl),
    acum(rhs.acum),
    comp(rhs.comp),
    beta_sqn(rhs.beta_sqn),
    m0(rhs.m0),
    pi_est(rhs.pi_est),
    dirich(rhs.dirich),
    epssum(rhs.epssum),
    sigmae(rhs.sigmae),
    sigmag(rhs.sigmag),
    mu(rhs.mu) {
    epsilon = (double*) _mm_malloc(size_t(im4*4) * sizeof(double), 64);
    check_malloc(epsilon, __LINE__, __FILE__);
    for (int i=0; i<im4*4; i++)
        epsilon[i] = rhs.epsilon[i];
    mave = (double*) _mm_malloc(size_t(M) * sizeof(double), 64);
    check_malloc(mave, __LINE__, __FILE__);
    msig = (double*) _mm_malloc(size_t(M) * sizeof(double), 64);
    check_malloc(msig, __LINE__, __FILE__);
    for (int i=0; i<M; i++) {
        mave[i] = rhs.mave[i];
        msig[i] = rhs.msig[i];
    }
    cass = (int*) _mm_malloc(size_t(K * G) * sizeof(int), 64);
    check_malloc(cass, __LINE__, __FILE__);
}

void Phenotype::print_cass() {
    printf("INFO   : cass for phenotype [...]\n");
    for (int i=0; i<G; i++) {
        printf("         %2d : ", i);
        for (int j=0; j<K; j++) {
            printf("%7d", cass[i*K+j]);
        }
        printf("\n");
    }
}

void Phenotype::print_cass(const std::vector<int>& mtotgrp) {
    printf("INFO   : cass for phenotype [...]\n");
    for (int i=0; i<G; i++) {
        printf("         group %2d: %7d | cass: ", i, mtotgrp.at(i));
        for (int j=0; j<K; j++) {
            printf("%7d", cass[i*K + j]);
        }
        printf("\n");
    }
}

void Phenotype::update_pi_est_dirichlet(const int g) {
    std::vector<double> tmp;
    double sum = 0.0;
    for (int i=0; i<K; i++) {
        double val = dist.rgamma((double)cass[g * K + i] + dirich[i], 1.0);
        set_pi_est(g, i, val);
        sum += val;
    }
    for (int i=0; i<K; i++)
        set_pi_est(g, i, get_pi_est(g, i) / sum);
}

double Phenotype::epsilon_sum() {
    double sum = 0.0;
    for (int i=0; i<N; i++) {
        sum += epsilon[i];
    }
    return sum;
}

double Phenotype::epsilon_sumsqr() {
    double sumsqr = 0.0;
    for (int i=0; i<N; i++) {
        sumsqr += epsilon[i] * epsilon[i];
    }
    return sumsqr;
}

void Phenotype::increment_beta_sqn(const int group, const double val) {
    beta_sqn.at(group) += val;
}

int Phenotype::get_marker_local_index(const int shuff_idx) {
    return midx[shuff_idx];
}

double Phenotype::sample_inv_scaled_chisq_rng(const double a, const double b) {
    return dist.inv_scaled_chisq_rng(a, b);
}

double Phenotype::sample_norm_rng(const double a, const double b) {
    return dist.norm_rng(a, b);
}

double Phenotype::sample_norm_rng() {
    //printf("sampling mu with epssum = %20.15f and sigmae = %20.15f; nonas = %d\n", epssum, sigmae, nonas);
    return dist.norm_rng(epssum / double(nonas), sigmae / double(nonas));
}

double Phenotype::sample_beta_rng(const double a, const double b) {
    return dist.beta_rng(a, b);
}

double Phenotype::sample_unif_rng() {
    return dist.unif_rng();
}

//void Phenotype::sample_sigmag_beta_rng() {
//    sigmag = dist.beta_rng(epssum / double(nonas), sigmae / double(nonas));
//}

void Phenotype::set_rng(const unsigned int seed) {
    dist.set_rng(seed);
}

void Phenotype::set_midx() {
    midx.clear();
    for (int i=0; i<M; ++i)
        midx.push_back(i);
}

void Phenotype::shuffle_midx() {
    boost::uniform_int<> unii(0, M-1);
    boost::variate_generator< boost::mt19937&, boost::uniform_int<> > generator(dist.get_rng(), unii);
    boost::range::random_shuffle(midx, generator);
}

// Add contributions to base epsilon
void Phenotype::update_epsilon(const double* dbeta, const unsigned char* bed) {

    const double bs_ = dbeta[0] * dbeta[2]; // lambda = dbeta / marker_sig
    //printf(" Phenotype::update_epsilon with dbeta = %20.15f, ave = %20.15f, bet/sig = %20.15f\n", dbeta[0], dbeta[1], bs_);

    __m256d mu = _mm256_set1_pd(-1.0 * dbeta[1]);
    __m256d bs = _mm256_set1_pd(bs_);
    __m256d eps4, deps4, lutna, luta, lutb;

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif                              

    for (int i=0; i<im4; i++) {
        eps4  = _mm256_load_pd(&epsilon[i*4]);
        lutna = _mm256_load_pd(&na_lut[mask4[i] * 4]);
        luta  = _mm256_load_pd(&dotp_lut_a[bed[i] * 4]);
        lutb  = _mm256_load_pd(&dotp_lut_b[bed[i] * 4]);
        deps4 = _mm256_fmadd_pd(mu, lutb, luta);
        deps4 = _mm256_mul_pd(deps4, bs);
        deps4 = _mm256_mul_pd(deps4, lutna);
        eps4  = _mm256_add_pd(eps4, deps4);
        _mm256_store_pd(&epsilon[i*4], eps4);
    }
}

// Assume all operations to be masked, so don't care about
// potential extra individuals from last byte of bed
void Phenotype::offset_epsilon(const double offset) {
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
    for (int i=0; i<im4*4; i++)
        epsilon[i] += offset;
}

void Phenotype::update_epsilon_sum() {
    __m256d eps4, sig4, lutna;
    __m256d sums = _mm256_set1_pd(0.0);

#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(addpd4:sums)
#endif
    for (int i=0; i<im4; i++) {
        eps4  = _mm256_load_pd(&epsilon[i*4]);
        lutna = _mm256_load_pd(&na_lut[mask4[i] * 4]);
        eps4  = _mm256_mul_pd(eps4, lutna);
        sums  = _mm256_add_pd(sums, eps4);
    }
    epssum = sums[0] + sums[1] + sums[2] + sums[3];
}

// Only depends on NAs 
void Phenotype::update_epsilon_sigma() {
    __m256d eps4, sig4, lutna;
    __m256d sume = _mm256_set1_pd(0.0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(addpd4:sume)
#endif
    for (int i=0; i<im4; i++) {
        eps4  = _mm256_load_pd(&epsilon[i*4]);
        lutna = _mm256_load_pd(&na_lut[mask4[i] * 4]);
        eps4  = _mm256_mul_pd(eps4, lutna);
        sig4  = _mm256_mul_pd(eps4, eps4);
        sume  = _mm256_add_pd(sume, sig4);
    }
    sigmae = (sume[0] + sume[1] + sume[2] + sume[3]) / double(nonas) * 0.5;
}


// Compute mean and associated standard deviation for markers
// for each of the phenotypes (stats are NA dependent)
// ! one byte of bed  contains information for 4 individuals
// ! one byte of phen contains information for 8 individuals
void PhenMgr::compute_markers_statistics(const unsigned char* bed, const int N, const int M, const int mbytes) {

    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for (auto& phen : get_phens()) {

        if (rank == 0)
            phen.print_info();

        const std::vector<unsigned char> mask4 = phen.get_mask4();

        double* mave = phen.get_mave();
        double* msig = phen.get_msig();

        //double start = MPI_Wtime();
#ifdef _OPENMP
#pragma omp parallel for
#endif 
        for (int i=0; i<M; i++) {
            __m256d suma = _mm256_set1_pd(0.0);
            __m256d sumb = _mm256_set1_pd(0.0);
            __m256d luta, lutb, lutna;
            size_t bedix = size_t(i) * size_t(mbytes);
            const unsigned char* bedm = &bed[bedix];
            for (int j=0; j<mbytes; j++) {
                luta  = _mm256_load_pd(&dotp_lut_a[bedm[j] * 4]);
                lutb  = _mm256_load_pd(&dotp_lut_b[bedm[j] * 4]);
                lutna = _mm256_load_pd(&na_lut[mask4[j] * 4]);
                luta  = _mm256_mul_pd(luta, lutna);
                lutb  = _mm256_mul_pd(lutb, lutna);
                suma  = _mm256_add_pd(suma, luta);
                sumb  = _mm256_add_pd(sumb, lutb);
            }
            double asum = suma[0] + suma[1] + suma[2] + suma[3];
            double bsum = sumb[0] + sumb[1] + sumb[2] + sumb[3];
            double avg  = round_dp(asum / bsum);

            __m256d vave = _mm256_set1_pd(-avg);
            __m256d sums = _mm256_set1_pd(0.0);
            for (int j=0; j<mbytes; j++) {
                luta  = _mm256_load_pd(&dotp_lut_a[bedm[j] * 4]);
                lutb  = _mm256_load_pd(&dotp_lut_b[bedm[j] * 4]);
                lutna = _mm256_load_pd(&na_lut[mask4[j] * 4]);
                luta  = _mm256_add_pd(luta, vave);    // - mu
                luta  = _mm256_mul_pd(luta, lutb);    // M -> 0.0
                luta  = _mm256_mul_pd(luta, lutna);   // NAs
                luta  = _mm256_mul_pd(luta, luta);    // ^2
                sums  = _mm256_add_pd(sums, luta);    // sum
            }
            double sig = round_dp(1.0 / sqrt((sums[0] + sums[1] + sums[2] + sums[3]) / (double(phen.get_nonas()) - 1.0)));
            mave[i] = avg;
            msig[i] = sig;
            //if (i<10)
            //    printf("marker %d: %20.15f +/- %20.15f, %20.15f / %20.15f\n", i, mave[i], msig[i], asum, bsum);
        }
        //double end = MPI_Wtime();
        //std::cout << "statistics took " << end - start << " seconds to run." << std::endl;
    }
}

void PhenMgr::display_markers_statistics(const int n) {
    for (auto& phen : get_phens()) {
        phen.print_info();
        for (int i=0; i<n; i++) {
            printf("avg for marker %d = %20.15f +/- %20.15f  1/sig = %20.15f\n", i, phen.get_mave()[i], phen.get_msig()[i], 1.0 / phen.get_msig()[i]);
        }
    }
}

void PhenMgr::print_info() {
    for (auto& phen : get_phens())
        phen.print_info();
}

void PhenMgr::read_phen_files(const Options& opt, const int N, const int M) {
    std::vector<std::string> phen_files = opt.get_phen_files();
    int pi = 0;
    for (auto fp = phen_files.begin(); fp != phen_files.end(); ++fp) {
        pi++;
        if (opt.verbosity_level(3))
            std::cout << "Reading phenotype file " << pi << ": " << *fp << std::endl;
        phens.emplace_back(*fp, opt, N, M);
        //std::cout << "emplace_back for pi " << pi << ": " << *fp << std::endl;
    }
}


// Read phenotype file assuming PLINK format:
// Family ID, Individual ID, Phenotype; One row per individual
void Phenotype::read_file(const Options& opt) {

    std::ifstream infile(filepath);
    std::string line;
    std::regex re("\\s+");

    double sum = 0.0;

    if (infile.is_open()) {
        int line_n = 0;
        nonas = 0, nas = 0;
        while (getline(infile, line)) {            
            int m4 = line_n % 4;
            if (m4 == 0)  mask4.push_back(0b00001111);

            std::sregex_token_iterator first{line.begin(), line.end(), re, -1}, last;
            std::vector<std::string> tokens{first, last};
            if (tokens[2] == "NA") {
                nas += 1;
                data.push_back(nan(0));
                if (opt.verbosity_level(3))
                    std::cout << " ... found NA on line " << line_n << ", m4 = " << m4 << " on byte " << int(line_n / 4) << std::endl;
                mask4.at(int(line_n / 4)) &= ~(0b1 << m4);
            } else {
                nonas += 1;
                data.push_back(atof(tokens[2].c_str()));
                sum += atof(tokens[2].c_str());
            }

            line_n += 1;

            if (opt.verbosity_level(3)) {
                if (line_n % 4 == 0 && line_n > 3 && line_n < 30) {
                    std::cout << "mask4[" << int(line_n / 4) - 1 << "] = " << unsigned(mask4.at(int(line_n / 4) - 1)) << std::endl;
                }
            }
        }
        infile.close();

        assert(nas + nonas == N);
        
        // Set last bits to 0 if ninds % 4 != 0
        const int m4 = line_n % 4;
        if (m4 != 0) {
            for (int i=m4; i<4; i++) {
                mask4.at(int(line_n / 4)) &= ~(0b1 << i);
            }
            //printf("line_n = %d\n", line_n);
            //printf("last byte starts for indiv %d\n", int(N/4)*4);
            //printf("set up to indiv %d\n", int(N/4 + 1) * 4);
            std::cout << "Setting last " << 4 - m4 << " bits to NAs" << std::endl;
            //std::cout << "fatal: missing implementation" << std::endl;
            //exit(1);
        }
        
        // Center and scale
        double avg = sum / double(nonas);
        if (opt.verbosity_level(3))
            printf("phen avg = %20.15f\n", avg);

        double sqn = 0.0;
        for (int i=0; i<data.size(); i++) {
            if (opt.verbosity_level(3) && i < 20)
                std::cout << data[i] - avg  << std::endl;
            if (! isnan(data[i])) {
                epsilon[i] = data[i] - avg;
                sqn += epsilon[i] * epsilon[i];
            } else {
                epsilon[i] = 0.0;
            }
        }
        sqn = sqrt(double(nonas-1) / sqn);
        if (opt.verbosity_level(3))
            printf("phen sqn = %20.15f\n", sqn);

        for (int i=0; i<data.size(); i++)
            epsilon[i] *= sqn;

    } else {
        std::cout << "FATAL: could not open phenotype file: " << filepath << std::endl;
        exit(EXIT_FAILURE);
    }
}

void Phenotype::print_info() const {
    printf("INFO   : %s has %d NAs and %d non-NAs.\n", get_filepath().c_str(), nas, nonas);
}
