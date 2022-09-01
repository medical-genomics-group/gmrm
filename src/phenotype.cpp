#pragma once
#include <iostream>
#include <fstream>
#include <regex>
#include <math.h>
#include <immintrin.h>
//#include <mpi.h>
//#include "utilities.hpp"
#include "phenotype.hpp"
#include "dotp_lut.hpp"
#include "na_lut.hpp"

#include <boost/range/algorithm.hpp>
#include <boost/random/uniform_int.hpp>
#define BOOST_UBLAS_NDEBUG 1
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/random.hpp>


#include <filesystem>

namespace fs = std::filesystem;


Phenotype::Phenotype(std::string fp, const Options& opt, const int N, const int M) :
    filepath(fp),
    N(N),
    M(M),
    im4(N % 4 == 0 ? N / 4 : N / 4 + 1),
    K(opt.get_nmixtures()),
    G(opt.get_ngroups()) {

    mave = (double*)_mm_malloc(size_t(M) * sizeof(double), 32);
    //check_malloc(mave, __LINE__, __FILE__);
    msig = (double*)_mm_malloc(size_t(M) * sizeof(double), 32);
    //check_malloc(msig, __LINE__, __FILE__);
    epsilon_ = (double*)_mm_malloc(size_t(im4 * 4) * sizeof(double), 32);
    //check_malloc(epsilon_, __LINE__, __FILE__);

    if (!opt.predict()) {
        //cass = (int*) _mm_malloc(g * k * sizeof(int), 32);
        //check_malloc(cass, __line__, __file__);

    //    betas.assign(m, 0.0);
    //    acum.assign(m, 0.0);
    //    muk.assign(k, 0.0);
    //    denom.resize(k - 1);
    //    logl.resize(k);
    //    comp.assign(m, 0);
    //    beta_sqn.resize(g);
    //    m0.resize(g);
    //    sigmag.resize(g);

    //    dirich.clear();
    //    for (int i=0; i<K; i++) dirich.push_back(1.0);
    //} else {
    //    set_prediction_filenames(opt.get_out_dir());
    //}

        //set_output_filenames(opt.get_out_dir());
        read_file(opt); // Files are read, epsilons, nas, nonas, data are created
        cout << "Phenotype file read." << "\n";
    }
}

//int Phenotype::get_m0_sum() {
//    int sum = 0;
//    for (int i=0; i<m0.size(); i++)
//        sum += m0.at(i);
//    return sum;
//}

// Set input filenames based on input phen file (as per output)
//void Phenotype::set_prediction_filenames(const std::string out_dir) {
//    fs::path pphen = filepath;
//    fs::path base  = out_dir;
//    base /= pphen.stem();
//    fs::path pibet = base;
//    pibet.replace_extension(".bet");
//    inbet_fp = pibet.string();
//    //std::cout << "inbet_fp = " << inbet_fp << std::endl;
//
//    fs::path pmlma = base;
//    pmlma += ".mlma";
//    outmlma_fp = pmlma.string();
//}

//void Phenotype::set_output_filenames(const std::string out_dir) {
//    fs::path pphen = filepath;
//    fs::path base  = out_dir;
//    base /= pphen.stem();
//    fs::path pbet = base;
//    pbet += ".bet";
//    fs::path pcpn = base;
//    pcpn += ".cpn";
//    fs::path pcsv = base;
//    pcsv += ".csv";
//
//    outbet_fp = pbet.string();
//    outcpn_fp = pcpn.string();
//    outcsv_fp = pcsv.string();
//}

// Input and output
//void Phenotype::open_prediction_files() {
//
//    check_mpi(MPI_File_open(MPI_COMM_WORLD,
//                            get_inbet_fp().c_str(),
//                            MPI_MODE_RDONLY,
//                            MPI_INFO_NULL,
//                            get_inbet_fh()),
//              __LINE__, __FILE__);
//
//    check_mpi(MPI_File_open(MPI_COMM_WORLD,
//                            get_outmlma_fp().c_str(),
//                            MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL,
//                            MPI_INFO_NULL,
//                            get_outmlma_fh()),
//              __LINE__, __FILE__);
//}
//
//void Phenotype::close_prediction_files() {
//    check_mpi(MPI_File_close(get_inbet_fh()), __LINE__, __FILE__);
//    check_mpi(MPI_File_close(get_outmlma_fh()), __LINE__, __FILE__);
//}
//
//void Phenotype::delete_output_prediction_files() {
//    MPI_File_delete(get_outmlma_fp().c_str(), MPI_INFO_NULL);
//}
//
//void Phenotype::open_output_files() {
//    check_mpi(MPI_File_open(MPI_COMM_WORLD,
//                            get_outbet_fp().c_str(),
//                            MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL,
//                            MPI_INFO_NULL,
//                            get_outbet_fh()),
//              __LINE__, __FILE__);
//    check_mpi(MPI_File_open(MPI_COMM_WORLD,
//                            get_outcpn_fp().c_str(),
//                            MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL,
//                            MPI_INFO_NULL,
//                            get_outcpn_fh()),
//              __LINE__, __FILE__);
//    check_mpi(MPI_File_open(MPI_COMM_WORLD,
//                            get_outcsv_fp().c_str(),
//                            MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL,
//                            MPI_INFO_NULL,
//                            get_outcsv_fh()),
//              __LINE__, __FILE__);
//}
//
//void Phenotype::close_output_files() {
//    check_mpi(MPI_File_close(get_outbet_fh()), __LINE__, __FILE__);
//    check_mpi(MPI_File_close(get_outcpn_fh()), __LINE__, __FILE__);
//    check_mpi(MPI_File_close(get_outcsv_fh()), __LINE__, __FILE__);
//}
//
//void Phenotype::delete_output_files() {
//    MPI_File_delete(get_outbet_fp().c_str(), MPI_INFO_NULL);
//    MPI_File_delete(get_outcpn_fp().c_str(), MPI_INFO_NULL);
//    MPI_File_delete(get_outcsv_fp().c_str(), MPI_INFO_NULL);
//}
//
//void Phenotype::print_cass() {
//    printf("INFO   : cass for phenotype [...]\n");
//    for (int i=0; i<G; i++) {
//        printf("         %2d : ", i);
//        for (int j=0; j<K; j++) {
//            printf("%7d", cass[i*K+j]);
//        }
//        printf("\n");
//    }
//}
//
//void Phenotype::print_cass(const std::vector<int>& mtotgrp) {
//    printf("INFO   : cass for phenotype [...]\n");
//    for (int i=0; i<G; i++) {
//        printf("         group %2d: %7d | cass: ", i, mtotgrp.at(i));
//        for (int j=0; j<K; j++) {
//            printf("%7d", cass[i*K + j]);
//        }
//        printf("\n");
//    }
//}
//
//void Phenotype::update_pi_est_dirichlet(const int g) {
//    std::vector<double> tmp;
//    double sum = 0.0;
//    for (int i=0; i<K; i++) {
//        double val = dist_d.rgamma((double)cass[g * K + i] + dirich[i], 1.0);
//        set_pi_est(g, i, val);
//        sum += val;
//    }
//    for (int i=0; i<K; i++)
//        set_pi_est(g, i, get_pi_est(g, i) / sum);
//}

double Phenotype::epsilon_sum() {
    double* epsilon = get_epsilon();
    double sum = 0.0;
//#ifdef _OPENMP
//#pragma omp parallel for simd aligned(epsilon:32) reduction(+:sum)
//#endif
    for (int i=0; i<N; i++) {
        sum += epsilon[i];
    }
    return sum;
}

double Phenotype::epsilon_sumsqr() {
    double* epsilon = get_epsilon();
    double sumsqr = 0.0;
#ifdef _OPENMP
#pragma omp parallel for simd aligned(epsilon:32) reduction(+:sumsqr)
#endif
    for (int i=0; i<N; i++) {
        sumsqr += epsilon[i] * epsilon[i];
    }
    return sumsqr;
}

//void Phenotype::increment_beta_sqn(const int group, const double val) {
//    beta_sqn.at(group) += val;
//}
//
//int Phenotype::get_marker_local_index(const int shuff_idx) {
//    return midx[shuff_idx];
//}
//
//double Phenotype::sample_inv_scaled_chisq_rng(const double a, const double b) {
//    return dist_d.inv_scaled_chisq_rng(a, b);
//}
//
//double Phenotype::sample_norm_rng(const double a, const double b) {
//    return dist_d.norm_rng(a, b);
//}
//
//double Phenotype::sample_norm_rng() {
//    //printf("sampling mu with epssum = %20.15f and sigmae = %20.15f; nonas = %d\n", epssum, sigmae, nonas);
//    return dist_d.norm_rng(epssum / double(nonas), get_sigmae() / double(nonas));
//}
//
//double Phenotype::sample_beta_rng(const double a, const double b) {
//    return dist_d.beta_rng(a, b);
//}
//
//double Phenotype::sample_unif_rng() {
//    return dist_d.unif_rng();
//}
//
//void Phenotype::sample_for_free(const int n) {
//    for (int i=0; i<n; i++)
//        double fake = dist_d.unif_rng();
//}

//void Phenotype::sample_sigmag_beta_rng() {
//    sigmag = dist.beta_rng(epssum / double(nonas), sigmae / double(nonas));
//}

//void Phenotype::set_prng_m(const unsigned int seed) {
//    dist_m.set_prng(seed);
//}
//void Phenotype::set_prng_d(const unsigned int seed) {
//    dist_d.set_prng(seed);
//}
//
//void Phenotype::set_midx() {
//    midx.clear();
//    for (int i=0; i<M; ++i)
//        midx.push_back(i);
//}
//
//void Phenotype::shuffle_midx(const bool mimic_hydra) {
//    boost::uniform_int<> unii(0, M-1);
//    if (mimic_hydra) {
//        boost::variate_generator< boost::mt19937&, boost::uniform_int<> > generator(dist_d.get_rng(), unii);
//        boost::range::random_shuffle(midx, generator);
//    } else {
//        boost::variate_generator< boost::mt19937&, boost::uniform_int<> > generator(dist_m.get_rng(), unii);
//        boost::range::random_shuffle(midx, generator);
//    }
//}

// Update epsilon for concrete mrker of concrete phenotype
void Phenotype::update_epsilon(const double* dbeta, const unsigned char* bed) {

    const double bs_ = dbeta[0] * dbeta[2]; // lambda = dbeta / marker_sig
    const double mdb = dbeta[1];
    double* epsilon = get_epsilon();

    //#ifdef MANVECT
    //    const __m256d mu = _mm256_set1_pd(-1.0 * dbeta[1]);
    //    const __m256d bs = _mm256_set1_pd(bs_);
    //#ifdef __INTEL_COMPILER
    //    __assume_aligned(epsilon, 32);
    //    __assume_aligned(&na_lut, 32);
    //    __assume_aligned(&dotp_lut_a, 32);
    //    __assume_aligned(&dotp_lut_b, 32);
    //#endif
    //#ifdef _OPENMP
    //#pragma omp parallel for schedule(static)
    //#endif
    //    for (int i=0; i<im4; i++) {
    //        __m256d luta  = _mm256_load_pd(&dotp_lut_ab[bed[i] * 8]);
    //        __m256d lutb  = _mm256_load_pd(&dotp_lut_ab[bed[i] * 8 + 4]);
    //        __m256d tmp1  = _mm256_fmadd_pd(mu, lutb, luta);
    //        __m256d lutna = _mm256_load_pd(&na_lut[mask4[i] * 4]);
    //        __m256d tmp2  = _mm256_mul_pd(bs, tmp1);
    //        __m256d eps4  = _mm256_load_pd(&epsilon[i*4]);
    //        __m256d tmp3  = _mm256_fmadd_pd(lutna, tmp2, eps4);
    //        _mm256_store_pd(&epsilon[i*4], tmp3);
    //    }
    /*
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
    */

    //#else
    //
    //#ifdef _OPENMP
    //#pragma omp parallel for
    //#endif
    for (int i = 0; i < im4; i++) {
        const int bedi = bed[i] * 4;
        const int masi = mask4[i] * 4;
        //#ifdef _OPENMP
        //#pragma omp simd aligned(epsilon, dotp_lut_a, dotp_lut_b, na_lut : 32) simdlen(4)
        //#endif
        for (int j = 0; j < 4; j++) {
            double a = dotp_lut_a[bedi + j];
            double b = dotp_lut_b[bedi + j];
            double m = na_lut[masi + j];
            //std::cout << "eps was " << epsilon[i * 4 + j] << "\n";
            //std::cout << "added val " << (a - mdb * b) * bs_ * m << "\n";
            epsilon[i * 4 + j] += (a - mdb * b) * bs_ * m;
            //std::cout << "eps new " << epsilon[i * 4 + j] << "\n";
        }
    }

    //#endif
}



//// Add contributions to base epsilon
//void Phenotype::update_epsilon(const double* dbeta, const unsigned char* bed) {
//
//    const double bs_ = dbeta[0] * dbeta[2]; // lambda = dbeta / marker_sig
//    const double mdb = -dbeta[1];
//    //printf(" Phenotype::update_epsilon with dbeta = %20.15f, ave = %20.15f, bet/sig = %20.15f\n", dbeta[0], dbeta[1], bs_);
//
//    double* epsilon = get_epsilon();
//
////#ifdef MANVECT
////    const __m256d mu = _mm256_set1_pd(-1.0 * dbeta[1]);
////    const __m256d bs = _mm256_set1_pd(bs_);
////#ifdef __INTEL_COMPILER
////    __assume_aligned(epsilon, 32);
////    __assume_aligned(&na_lut, 32);
////    __assume_aligned(&dotp_lut_a, 32);
////    __assume_aligned(&dotp_lut_b, 32);
////#endif
////#ifdef _OPENMP
////#pragma omp parallel for schedule(static)
////#endif
////    for (int i=0; i<im4; i++) {
////        __m256d luta  = _mm256_load_pd(&dotp_lut_ab[bed[i] * 8]);
////        __m256d lutb  = _mm256_load_pd(&dotp_lut_ab[bed[i] * 8 + 4]);
////        __m256d tmp1  = _mm256_fmadd_pd(mu, lutb, luta);
////        __m256d lutna = _mm256_load_pd(&na_lut[mask4[i] * 4]);
////        __m256d tmp2  = _mm256_mul_pd(bs, tmp1);
////        __m256d eps4  = _mm256_load_pd(&epsilon[i*4]);
////        __m256d tmp3  = _mm256_fmadd_pd(lutna, tmp2, eps4);
////        _mm256_store_pd(&epsilon[i*4], tmp3);
////    }
///*
//#ifdef _OPENMP
//#pragma omp parallel for schedule(static)
//#endif
//    for (int i=0; i<im4; i++) {
//        eps4  = _mm256_load_pd(&epsilon[i*4]);
//        lutna = _mm256_load_pd(&na_lut[mask4[i] * 4]);
//        luta  = _mm256_load_pd(&dotp_lut_a[bed[i] * 4]);
//        lutb  = _mm256_load_pd(&dotp_lut_b[bed[i] * 4]);
//        deps4 = _mm256_fmadd_pd(mu, lutb, luta);
//        deps4 = _mm256_mul_pd(deps4, bs);
//        deps4 = _mm256_mul_pd(deps4, lutna);
//        eps4  = _mm256_add_pd(eps4, deps4);
//        _mm256_store_pd(&epsilon[i*4], eps4);
//    }
//*/
//
////#else
////
////#ifdef _OPENMP
////#pragma omp parallel for
////#endif
//    for (int i=0; i<im4; i++) {
//        const int bedi = bed[i] * 4;
//        const int masi = mask4[i] * 4;
////#ifdef _OPENMP
////#pragma omp simd aligned(epsilon, dotp_lut_a, dotp_lut_b, na_lut : 32) simdlen(4)
////#endif
//        for (int j=0; j<4; j++) {
//            double a = dotp_lut_a[bedi + j];
//            double b = dotp_lut_b[bedi + j];
//            double m = na_lut[masi + j];
//            epsilon[i*4 + j] += (mdb * b + a) * bs_ * m;
//        }
//    }
//
////#endif
//}

void Phenotype::offset_epsilon(const double offset) {

    double* epsilon = get_epsilon();

//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
    for (int i=0; i<im4; i++) {
        const int masi = mask4[i] * 4;
//#ifdef _OPENMP
//#pragma omp simd aligned(epsilon, na_lut : 32) simdlen(4)
//#endif
        for (int j = 0; j < 4; j++) {
            epsilon[i * 4 + j] += offset * na_lut[masi + j];

        };
    };
}

/*
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
*/

// Only depends on NAs
double Phenotype::update_epsilon_sigma() {
    double* epsilon = get_epsilon();
#ifdef MANVECT
    __m256d sume = _mm256_set1_pd(0.0);
#ifdef _OPENMP
#pragma omp parallel for schedule(static) reduction(addpd4:sume)
#endif
    for (int i = 0; i < im4; i++) {
        __m256d eps4 = _mm256_load_pd(&epsilon[i * 4]);
        __m256d lutna = _mm256_load_pd(&na_lut[mask4[i] * 4]);
        eps4 = _mm256_mul_pd(eps4, lutna);
        eps4 = _mm256_mul_pd(eps4, eps4);
        sume = _mm256_add_pd(sume, eps4);
    }
    set_sigmae((sume[0] + sume[1] + sume[2] + sume[3]) / double(nonas) * 0.5);
#else
    double sigmae = 0.0;
#ifdef _OPENMP
#pragma omp parallel for simd reduction(+:sigmae)
#endif
    for (int i = 0; i < im4; i++) {
        for (int j = 0; j < 4; j++) {
            sigmae += epsilon[i * 4 + j] * epsilon[i * 4 + j] * na_lut[mask4[i] * 4 + j];
        };
    };
    //set_sigmae(sigmae / double(nonas) * 0.5);

    return sigmae / double(nonas) * 0.5;

#endif
}


// Compute mean and associated standard deviation for markers
// for each of the phenotypes (stats are NA dependent)
// ! one byte of bed  contains information for 4 individuals
// ! one byte of phen contains information for 8 individuals
void PhenMgr::compute_markers_statistics(const unsigned char* bed, const int N, const int M, const int mbytes) {

    //int rank;
    //MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for (auto& phen : get_phens()) {

        //if (rank == 0)
        //    phen.print_info();

        const std::vector<unsigned char> mask4 = phen.get_mask4();
        const int im4 = phen.get_im4();

        double* mave = phen.get_mave();
        double* msig = phen.get_msig();

        //double start = MPI_Wtime();
//#ifdef MANVECT
//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
//        for (int i=0; i<M; i++) {
//            __m256d suma = _mm256_set1_pd(0.0);
//            __m256d sumb = _mm256_set1_pd(0.0);
//            __m256d luta, lutb, lutna;
//            size_t bedix = size_t(i) * size_t(mbytes);
//            const unsigned char* bedm = &bed[bedix];
//            for (int j=0; j<mbytes; j++) {
//                luta  = _mm256_load_pd(&dotp_lut_a[bedm[j] * 4]);
//                lutb  = _mm256_load_pd(&dotp_lut_b[bedm[j] * 4]);
//                lutna = _mm256_load_pd(&na_lut[mask4[j] * 4]);
//                luta  = _mm256_mul_pd(luta, lutna);
//                lutb  = _mm256_mul_pd(lutb, lutna);
//                suma  = _mm256_add_pd(suma, luta);
//                sumb  = _mm256_add_pd(sumb, lutb);
//            }
//            double asum = suma[0] + suma[1] + suma[2] + suma[3];
//            double bsum = sumb[0] + sumb[1] + sumb[2] + sumb[3];
//            double avg  = asum / bsum;
//
//            __m256d vave = _mm256_set1_pd(-avg);
//            __m256d sums = _mm256_set1_pd(0.0);
//            for (int j=0; j<mbytes; j++) {
//                luta  = _mm256_load_pd(&dotp_lut_a[bedm[j] * 4]);
//                lutb  = _mm256_load_pd(&dotp_lut_b[bedm[j] * 4]);
//                lutna = _mm256_load_pd(&na_lut[mask4[j] * 4]);
//                luta  = _mm256_add_pd(luta, vave);    // - mu
//                luta  = _mm256_mul_pd(luta, lutb);    // M -> 0.0
//                luta  = _mm256_mul_pd(luta, lutna);   // NAs
//                luta  = _mm256_mul_pd(luta, luta);    // ^2
//                sums  = _mm256_add_pd(sums, luta);    // sum
//            }
//            double sig = 1.0 / sqrt((sums[0] + sums[1] + sums[2] + sums[3]) / (double(phen.get_nonas()) - 1.0));
//            mave[i] = avg;
//            msig[i] = sig;
//            //if (i<10)
//            //    printf("marker %d: %20.15f +/- %20.15f, %20.15f / %20.15f\n", i, mave[i], msig[i], asum, bsum);
//        }
//#else
//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
        for (int i=0; i<M; i++) {
            size_t bedix = size_t(i) * size_t(mbytes);
            const unsigned char* bedm = &bed[bedix];
            double suma = 0.0;
            double sumb = 0.0;
            for (int j=0; j<im4; j++) {
                for (int k=0; k<4; k++) {
                    suma += dotp_lut_a[bedm[j] * 4 + k] * na_lut[mask4[j] * 4 + k];
                    sumb += dotp_lut_b[bedm[j] * 4 + k] * na_lut[mask4[j] * 4 + k];
                }
            }
            mave[i] = suma / sumb;
            double sumsqr = 0.0;
            for (int j=0; j<im4; j++) {
                for (int k=0; k<4; k++) {
                    double val = (dotp_lut_a[bedm[j] * 4 + k] - mave[i]) * dotp_lut_b[bedm[j] * 4 + k] * na_lut[mask4[j] * 4 + k];
                    sumsqr += val * val;
                }
            }
            msig[i] = 1.0 / sqrt(sumsqr / (double(phen.get_nonas()) - 1.0));
        }

//#endif
        //double end = MPI_Wtime();
        //std::cout << "statistics took " << end - start << " seconds to run." << std::endl;
    }
}

void PhenMgr::display_markers_statistics(const int n) {
    for (auto& phen : get_phens()) {
        for (int i=0; i<n; i++) {
            printf("avg for marker %d = %20.15f +/- %20.15f  1/sig = %20.15f\n", i, phen.get_mave()[i], phen.get_msig()[i], 1.0 / phen.get_msig()[i]);
        }
    }
}

void PhenMgr::read_phen_files(const Options& opt, const int N1, const int M1) {
    N = N1;
    M = M1;
    G = opt.get_ngroups();
    K = opt.get_nmixtures();
    std::vector<std::string> phen_files = opt.get_phen_files();
    int pi = 0;
    for (auto fp = phen_files.begin(); fp != phen_files.end(); ++fp) {
        pi++;
        if (opt.verbosity_level(3))
            std::cout << "Reading phenotype file " << pi << ": " << *fp << std::endl;
        phens.emplace_back(*fp, opt, N1, M1); // A new phenotype is added and constructed at the same time
    }
}


// Read phenotype file assuming PLINK format:
// Family ID, Individual ID, Phenotype; One row per individual
void Phenotype::read_file(const Options& opt) {

    std::ifstream infile(filepath);
    std::string line;
    std::regex re("\\s+");

    double* epsilon = get_epsilon();
    double sum = 0.0;

    if (infile.is_open()) {
        int line_n = 0;
        nonas = 0, nas = 0;
        while (getline(infile, line)) {
            //std::cout << line << std::endl;
            int m4 = line_n % 4;
            if (m4 == 0)  mask4.push_back(0b00001111);

            std::sregex_token_iterator first{line.begin(), line.end(), re, -1}, last;
            std::vector<std::string> tokens{first, last};

            if (tokens[2] == "NA") {
                nas += 1;
                data.push_back(std::numeric_limits<double>::max());
                if (opt.verbosity_level(2)) {
                    std::cout << " ... found NA on line " << line_n << ", m4 = " << m4 << " on byte " << int(line_n / 4) << std::endl;
                    fflush(stdout);
                }
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
        //printf("phen avg = %20.15f\n", avg);

        double sqn = 0.0;
        for (int i=0; i<data.size(); i++) {
            if (opt.verbosity_level(3) && i < 10)
                std::cout << data[i] - avg  << std::endl;
            if (data[i] == std::numeric_limits<double>::max()) {
                epsilon[i] = 0.0;
            } else {
                epsilon[i] = data[i] - avg;
                sqn += epsilon[i] * epsilon[i];
            }
        }
        sqn = sqrt(double(nonas-1) / sqn);
        if (opt.verbosity_level(3))
            printf("phen sqn = %20.15f\n", sqn);

        for (int i = 0; i < data.size(); i++) {
            epsilon[i] *= sqn;
            //std::cout << epsilon[i] << "\n";
        };
    } else {
        std::cout << "FATAL: could not open phenotype file: " << filepath << std::endl;
        exit(EXIT_FAILURE);
    }
}


void Phenotype::set_nas_to_zero(double* y, const int N) {
    assert(N % 4 == 0);
    for (int j=0; j<N/4; j++) {
        for (int k=0; k<4; k++) {
            if (na_lut[mask4[j] * 4 + k] == 0.0) {
				assert(y[j*4 + k] == 0.0);
                //printf("found NA on %d. Correct? %20.15f\n", j * 4 + k, y[j*4 + k]);
                //y[j*4 + k] = 0.0;
            }
        }
    }
}



// ***********Multitriate***********


void PhenMgr::set_prng_m(const unsigned int seed) {
    dist_m.set_prng(seed);
}
void PhenMgr::set_prng_d(const unsigned int seed) {
    dist_d.set_prng(seed);
}

// Create initial identity covariance matrix for the groups
void PhenMgr::set_initial_sigmag(const int g, const int q, const int p) {

    for (int i = 0; i < g; i++)
    {
        std::vector<std::vector<double>> v2d;
        for (int j = 0; j < q; j++)
        {
            std::vector<double> v1d;
            for (int k = 0; k < q; k++)
            {
                if (j == k) {
                    v1d.push_back((double)1 / p);
                }
                else {
                    v1d.push_back(0.0);
                };
            };
            v2d.push_back(v1d);
        };
        sigmag.push_back(v2d);
    };
}

// Create initial vector for pi priors for all markers
void PhenMgr::set_pi_est(double val, const int p) {

    for (int i = 0; i < p; i++)
    {
        pi_est.push_back(val);
    };
}


// Create initial covariance matrix for the residuals
void PhenMgr::set_initial_sigmae(const int q, const int n) {
    for (int j = 0; j < q; j++)
    {
        std::vector<double> v1d;
        for (int k = 0; k < q; k++)
        {
            if (j == k) {
                v1d.push_back((double)1 / n);
            }
            else {
                v1d.push_back(0.0);
            };
        };
        sigmae.push_back(v1d);
    };
}

// Create initial vector of means for the residuals
void PhenMgr::set_initial_epsilon_mu(const int q) {
    for (int j = 0; j < q; j++)
    {
        mu_epsilon.push_back(0.0);
    };
}

// Set newmeans for epsilon
void PhenMgr::set_epsilon_mu(std::vector<double> val) {
    for (int j = 0; j < val.size(); j++)
    {
        mu_epsilon.at(j) = val.at(j);
    };
}


// Get sample means for all epsilons
std::vector<double> PhenMgr::get_epsilon_means() {
    std::vector<double> means;
    for (auto& phen : get_phens()) {
        means.push_back(phen.epsilon_sum() / double(phen.get_nonas()));
    };
    return means;
}


// Sample MVN
std::vector<double> PhenMgr::sample_mvnorm_rng(std::vector<double> mean, std::vector<std::vector<double>> sigma, int q) {

    std::vector<double> a;
    for (int i = 0; i < sigma.size(); ++i) {
        for (int j = 0; j < sigma[0].size(); ++j) {
            a.push_back(sigma[j][i]);
        };
    };

    std::vector<double> mus;
    dist_d.mvnorm_rng(q, 1, a, mean, 2022, mus);

    return mus;

}

// Setting markers indices
void PhenMgr::set_midx() {
    midx.clear();
    for (int i = 0; i < M; ++i)
        midx.push_back(i);
}

// Shuffling markers
void PhenMgr::shuffle_midx(const bool mimic_hydra) {
    boost::uniform_int<> unii(0, M - 1);
    if (mimic_hydra) {
        boost::variate_generator< boost::mt19937&, boost::uniform_int<> > generator(dist_d.get_rng(), unii);
        boost::range::random_shuffle(midx, generator);
    }
    else {
        boost::variate_generator< boost::mt19937&, boost::uniform_int<> > generator(dist_m.get_rng(), unii);
        boost::range::random_shuffle(midx, generator);
    }
}

// Get local index of marker
int PhenMgr::get_marker_local_index(const int shuff_idx) {
    return midx[shuff_idx];
}

// Setting the vector to track number of times each marker has been processed
void PhenMgr::set_acum() {
    for (int i = 0; i < M; ++i)
        acum.push_back(0);
}

// Set initial matrix of bets of size pxq
void PhenMgr::set_initial_betas(const int q) {
    for (int j = 0; j < M; j++)
    {
        std::vector<double> v1d;
        for (int k = 0; k < q; k++)
        {
            v1d.push_back(0.0);
        };
        betas.push_back(v1d);
    };
}

// Creation of inverse of matrix with boost
std::vector<std::vector<double>> PhenMgr::InvertMatrix(std::vector<std::vector<double>> input, int q) {
    using namespace boost::numeric::ublas;
    typedef permutation_matrix<std::size_t> pmatrix;

    // create a working copy of the input
    matrix<double> A(q,q);
    for (size_t i = 0; i < A.size1(); i++)
    {
        for (size_t j = 0; j < A.size2(); j++)
        {
            A(i, j) = input.at(i).at(j);
        };
    };

    // create the inverse ublas matrix
    matrix<double> inverse(q,q);
    // create a permutation matrix for the LU-factorization
    pmatrix pm(A.size1());

    int res = lu_factorize(A, pm);
    
    // create identity matrix of "inverse"
    inverse.assign(identity_matrix<double>(A.size1()));

    // backsubstitute to get the inverse
    lu_substitute(A, pm, inverse);

    // create 2D vector to put inverse matrix elements
    std::vector<std::vector<double>> inv(q, std::vector<double>(q));
    for (size_t i = 0; i < inverse.size1(); i++)
    {
        for (size_t j = 0; j < inverse.size2(); j++)
        {
            inv.at(i).at(j) = inverse(i, j);
        };
    };

    return inv;
}

// Multiplication of matrix by scalar
std::vector<std::vector<double>> PhenMgr::multiply_matrix_scalar(std::vector<std::vector<double>> v, double k) {
    std::vector<std::vector<double>> mult(v.size(), std::vector<double>(v.size()));
    
    for (int i = 0; i < v.size(); ++i) {
        for (int j = 0; j < v[0].size(); ++j) {
            mult[i][j] = v[i][j] * k;
        };
    };
    return mult;
}

// Multiplication of vector by scalar
std::vector<double> PhenMgr::multiply_vector_scalar(std::vector<double> v, int k) {
    std::vector<double> mult(v.size());

    for (int i = 0; i < v.size(); ++i) {
        mult[i] = v[i] * k;
    };
    return mult;
}

// Sum of matrices
std::vector<std::vector<double>> PhenMgr::sum_matrices(std::vector<std::vector<double>> m1, std::vector<std::vector<double>> m2) {
    std::vector<std::vector<double>> sum(m1.size(), std::vector<double>(m1.size()));
    
    for (int i = 0; i < m1.size(); ++i) {
        for (int j = 0; j < m1[0].size(); ++j) {
            sum[i][j] = m1[i][j] + m2[i][j];
        };
    };
    return sum;
}

// Sum of vectors
std::vector<double> PhenMgr::sum_vectors(std::vector<double> v1, std::vector<double> v2) {
    std::vector<double> sum;

    for (int i = 0; i < v1.size(); i++) {
        sum.push_back(v1[i] + v2[i]);
    };
    return sum;
}

// Vector matrix product
std::vector<double> PhenMgr::vec_matrix_product(std::vector<double> v1, std::vector<std::vector<double>> m) {
    std::vector<double> v2(m[0].size());

    for (int j = 0; j <m[0].size(); ++j)
    {
        for (int k = 0; k < v1.size(); ++k)
        {
            v2[j] += v1[k] * m[k][j];
        };
    };

    return v2;
}


// Matrix vector product
std::vector<double> PhenMgr::matrix_vec_product(std::vector<std::vector<double>> m, std::vector<double> v1) {
    std::vector<double> v2(m.size());

    for (int j = 0; j < m.size(); ++j)
    {
        for (int k = 0; k < v1.size(); ++k)
        {
            v2[j] += m[j][k] * v1[k];
        };
    };

    return v2;
}

// Dot product of two vectors
double PhenMgr::vec_vec_product(std::vector<double> v1, std::vector<double> v2) {
    double sum = 0.0;

    for (int j = 0; j < v1.size(); j++)
    {
        sum += v1[j] * v2[j];
    };

    return sum;
}


double PhenMgr::determinant(std::vector<std::vector<double>> m, const int n) {
    double det = 0.0;
    std::vector<std::vector<double>> submatrix(n, std::vector<double>(n));
    if (n == 2)
        return ((m[0][0] * m[1][1]) - (m[1][0] * m[0][1]));
    else {
        for (int x = 0; x < n; x++) {
            int subi = 0;
            for (int i = 1; i < n; i++) {
                int subj = 0;
                for (int j = 0; j < n; j++) {
                    if (j == x)
                        continue;
                    submatrix[subi][subj] = m[i][j];
                    subj++;
                }
                subi++;
            }
            det = det + (pow(-1, x) * m[0][x] * determinant(submatrix, n - 1));
        }
    }
    return det;
}

// Matrix multiplication
std::vector<std::vector<double>> PhenMgr::mul_matrices(std::vector<std::vector<double>> m1, std::vector<std::vector<double>> m2) {
    std::vector<std::vector<double>> mat_mul(m1[0].size(), std::vector<double>(m1[0].size()));
    
    for (int i = 0; i < m1[0].size(); i++)
    {
        for (int j = 0; j < m2[0].size(); j++)
        {
            for (int k = 0; k < m1.size(); k++)
            {
                mat_mul[i][j] += m1[k][i] * m2[k][j];
            };
        };
    };

    return mat_mul;
}

// Sample from uniform
double PhenMgr::sample_unif_rng() {
    return dist_d.unif_rng();
}


// Sample from InverseWishart
std::vector<std::vector<double>> PhenMgr::sample_inv_wishart_rng(const int q, const int dof, std::vector<std::vector<double>> precision_matrix) {

    std::vector<std::vector<double>> wis(q, std::vector<double>(q));

    std::vector<std::vector<double>> precision_matrix_inv = InvertMatrix(precision_matrix, q);

    std::vector<double> a;
    for (int i = 0; i < precision_matrix.size(); ++i) {
        for (int j = 0; j < precision_matrix[0].size(); ++j) {
            a.push_back(precision_matrix_inv[j][i]);
        };
    };

    dist_d.inv_wishart_rng(q, dof, a, wis);

    return wis;
}

// Sample from Beta
double PhenMgr::sample_beta(const int a, const int b) {
    return dist_d.beta_rng(a, b);
};

// Set output files for covariances for residuals and effects
void PhenMgr::set_output_filenames(const std::string out_dir) {
    fs::path pphen_rescov = "/result_rescov.csv";
    fs::path base_rescov = out_dir;
    base_rescov += pphen_rescov;
    outcsv_rescov_fp = base_rescov.string();

    fs::path pphen_betascov = "/result_betascov.csv";
    fs::path base_betascov = out_dir;
    base_betascov += pphen_betascov;
    outcsv_betascov_fp = base_betascov.string();

    fs::path pphen_betas = "/result_betas.csv";
    fs::path base_betas = out_dir;
    base_betas += pphen_betas;
    outcsv_betas_fp = base_betas.string();
}