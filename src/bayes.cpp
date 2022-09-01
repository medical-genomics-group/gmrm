#pragma once
#include <iostream>
#include <fstream>
#include <iterator>
#include <limits.h>
#include <cmath>
#include <chrono>
#include <immintrin.h>
#include <omp.h>
#include "bayes.hpp"
//#include "utilities.hpp"
#include "dotp_lut.hpp"
#include "na_lut.hpp"
//#include "xfiles.hpp"
#include <boost/math/special_functions/gamma.hpp>


//void Bayes::predict() {
//
//    MPI_Barrier(MPI_COMM_WORLD);
//    double ts = MPI_Wtime();
//
//    check_openmp();
//
//    MPI_Status status;
//    MPI_Offset file_size = 0;
//    MPI_File*  fh;
//
//    cross_bim_files();
//
//    int pidx = 0;
//    for (auto& phen : pmgr.get_phens()) {
//        pidx += 1;
//
//        phen.delete_output_prediction_files();
//        phen.open_prediction_files();
//
//        const std::vector<unsigned char> mask4 = phen.get_mask4();
//        const int im4 = phen.get_im4();
//
//        fh = phen.get_inbet_fh();
//        check_mpi(MPI_File_get_size(*fh, &file_size), __LINE__, __FILE__);
//        //printf("file_size = %u B\n", file_size);
//
//        // First element of the .bet is the total number of processed markers
//        // Then: iteration (uint) beta (double) for all markers
//        uint Mtot_ = 0;
//        MPI_Offset betoff = size_t(0);
//        check_mpi(MPI_File_read_at_all(*fh, betoff, &Mtot_, 1, MPI_UNSIGNED, &status), __LINE__, __FILE__);
//        if (Mtot_ != m_refrsid.size()) {
//            printf("Mismatch between expected and Mtot read from .bet file: %lu vs %d\n", rsid.size(), Mtot_);
//            MPI_Abort(MPI_COMM_WORLD, 1);
//        }
//
//        assert((file_size - sizeof(uint)) % (Mtot_ * sizeof(double) + sizeof(uint)) == 0);
//        uint niter = (file_size - sizeof(uint)) / (Mtot_ * sizeof(double) + sizeof(uint));
//        if (rank == 0)
//            printf("INFO   : Number of recorded iterations in .bet file %d: %u\n", pidx-1, niter);
//
//        double* beta_sum = (double*) _mm_malloc(size_t(Mtot_) * sizeof(double), 32);
//        check_malloc(beta_sum, __LINE__, __FILE__);
//        for (int i=0; i<Mtot_; i++) beta_sum[i] = 0.0;
//
//        double* beta_it = (double*) _mm_malloc(size_t(Mtot_) * sizeof(double), 32);
//        check_malloc(beta_it, __LINE__, __FILE__);
//
//        uint start_iter = 0;
//        //EO: use this one to speed up testing (avoids to read the entire bet history)
//        //if (niter > 3) start_iter = niter - 3;
//
//        for (uint i=start_iter; i<niter; i++) {
//            betoff
//                = sizeof(uint) // Mtot
//                + (sizeof(uint) + size_t(Mtot_) * sizeof(double)) * size_t(i)
//                + sizeof(uint);
//            check_mpi(MPI_File_read_at_all(*fh, betoff, beta_it, Mtot_, MPI_DOUBLE, &status), __LINE__, __FILE__);
//            for (int j=0; j<Mtot_;j++)
//                beta_sum[j] += beta_it[j];
//        }
//
//        for (int j=0; j<Mtot_;j++)
//            beta_sum[j] /= double(niter);
//
//        fflush(stdout);
//        double t_1 = MPI_Wtime();
//        MPI_Barrier(MPI_COMM_WORLD);
//        //if (rank >= 0)
//        //    printf("INFO   : intermediate time 1 = %.2f seconds.\n", t_1 - ts);
//
//
//        double* g_k = (double*) _mm_malloc(size_t(im4*4) * sizeof(double), 32);
//        check_malloc(g_k, __LINE__, __FILE__);
//        for (int i=0; i<im4*4; i++) g_k[i] = 0.0;
//
//#ifdef _OPENMP
//#pragma omp parallel for schedule(static)
//#endif
//        for (int mrki=0; mrki<M; mrki++) {
//
//            // Get rsid from current
//            int mglo = S + mrki;
//            std::string id = rsid.at(mglo);
//
//            // Skip markers with no corresponding rsid in reference bim file
//            if (m_refrsid.find(id) == m_refrsid.end()) {
//                //printf("%d -> %d = %s not found in reference bim\n", mrki, mglo, id.c_str());
//                continue;
//            }
//
//            int rmglo = m_refrsid.find(id)->second; // global marker index in reference bed
//
//            size_t bedix = size_t(mrki) * size_t(mbytes);
//            const unsigned char* bedm = &bed_data[bedix];
//
//            double mave = phen.get_marker_ave(mrki);
//            double msig = phen.get_marker_sig(mrki);
//
//            for (int j=0; j<im4; j++) {
//                for (int k=0; k<4; k++) {
//                    double val = (dotp_lut_a[bedm[j] * 4 + k] - mave) * dotp_lut_b[bedm[j] * 4 + k] * na_lut[mask4[j] * 4 + k] * msig;
//#ifdef _OPENMP
//#pragma omp atomic update
//#endif
//                    g_k[j*4+k] += val * beta_sum[mglo];
//                }
//            }
//        }
//
//        fflush(stdout);
//        double t_2 = MPI_Wtime();
//        MPI_Barrier(MPI_COMM_WORLD);
//        //if (rank >= 0)
//        //    printf("INFO   : intermediate time 2 = %.2f seconds.\n", t_2 - t_1);
//
//        double* g = (double*) _mm_malloc(size_t(im4*4) * sizeof(double), 32);
//        check_malloc(g, __LINE__, __FILE__);
//
//        check_mpi(MPI_Allreduce(g_k, g, im4*4, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD), __LINE__, __FILE__);
//
//        double* y_k = (double*) _mm_malloc(size_t(im4*4) * sizeof(double), 32);
//        check_malloc(y_k, __LINE__, __FILE__);
//
//        phen.get_centered_and_scaled_y(y_k);
//
//        //EO: already done when computing original epsilon when reading .phen
//	//phen.set_nas_to_zero(y_k, im4*4);
//        // NAs are 0.0 in all, so should be safe
//        for (int i=0; i<im4*4; i++)
//            y_k[i] -= (g[i] - g_k[i]);
//
//        double sigma = 0.0;
//        for (int i=0; i<N; i++)
//            sigma += y_k[i] * y_k[i];
//        sigma /= phen.get_nonas();
//        //printf("### r: %d sigma = %20.15f\n", rank, sigma);
//
//        MPI_File* mlma_fh = phen.get_outmlma_fh();
//
//        double* Beta  = (double*) _mm_malloc(size_t(M) * sizeof(double), 32);
//        check_malloc(Beta, __LINE__, __FILE__);
//        double* Tdist = (double*) _mm_malloc(size_t(M) * sizeof(double), 32);
//        check_malloc(Tdist, __LINE__, __FILE__);
//        double* Se    = (double*) _mm_malloc(size_t(M) * sizeof(double), 32);
//        check_malloc(Se, __LINE__, __FILE__);
//        double* Pval  = (double*) _mm_malloc(size_t(M) * sizeof(double), 32);
//        check_malloc(Pval, __LINE__, __FILE__);
//
//
//#ifdef _OPENMP
//#pragma omp parallel for
//#endif
//        for (int mrki=0; mrki<M; mrki++) {
//
//            // Global index in current .bim
//            int mglo = S + mrki;
//            std::string id = rsid.at(mglo);
//
//            // Skip markers with no corresponding rsid in reference bim file
//            if (m_refrsid.find(id) == m_refrsid.end()) {
//                //printf("%d -> %d = %s not found in reference bim\n", mrki, mglo, id.c_str());
//                continue;
//            }
//
//            // Global marker index in reference bed
//            int rmglo = m_refrsid.find(id)->second;
//
//            size_t bedix = size_t(mrki) * size_t(mbytes);
//            const unsigned char* bedm = &bed_data[bedix];
//
//            double mave = phen.get_marker_ave(mrki);
//            double msig = phen.get_marker_sig(mrki);
//
//            double xtx = 0.0;
//            double xty = 0.0;
//            double chk = 0.0;
//            for (int j=0; j<im4; j++) {
//                for (int k=0; k<4; k++) {
//                    double val = dotp_lut_a[bedm[j] * 4 + k] * dotp_lut_b[bedm[j] * 4 + k] * na_lut[mask4[j] * 4 + k];
//                    xtx += val * val;
//                    xty += val * y_k[j*4 + k];
//                }
//            }
//
//            double beta  = xty / xtx;
//            double tdist = xty / sqrt(sigma * xtx);
//            double se    = beta / tdist;
//            double pval  = 1.0 - boost::math::gamma_p(0.5, tdist * tdist * 0.5);
//
//            Beta[mrki]  = beta;
//            Tdist[mrki] = tdist;
//            Se[mrki]    = se;
//            Pval[mrki]  = pval;
//
//            //if (mrki < 5) {
//            //    printf("r: %3d  m: %8d %8d (%5s) %8d  beta=%20.15f  se=%20.15f  tdist=%20.15f pval=%20.15f   xtx=%10.1f xty=%15.6f\n", rank, mrki, mglo, id.c_str(), rmglo, beta, se, tdist, pval, xtx, xty);
//            //}
//        }
//
//        fflush(stdout);
//        double t_3 = MPI_Wtime();
//        MPI_Barrier(MPI_COMM_WORLD);
//        //if (rank >= 0)
//        //    printf("INFO   : intermediate time 3 = %.2f seconds.\n", t_3 - t_2);
//
//        const int LLEN = 123 + 1;
//        char* todump  = (char*) _mm_malloc(size_t(LLEN) * size_t(M) * sizeof(char), 32);
//        check_malloc(todump, __LINE__, __FILE__);
//
//        int n_rem = 0;
//        for (int mrki=0; mrki<M; mrki++) {
//            int mglo = S + mrki;
//            std::string id = rsid.at(mglo);
//            if (m_refrsid.find(id) == m_refrsid.end()) {
//                //|| id.compare("rs12562034") == 0m|| id.compare("rs188466450") == 0)
//                printf("WARNING: marker id %s excluded -- no match\n", id.c_str());
//                n_rem++;
//                continue;
//            }
//            int rmglo = m_refrsid.find(id)->second;
//            int cx = snprintf(&todump[(mrki - n_rem) * (LLEN - 1)], LLEN, "%20s %8d %8d %20.15f %20.15f %20.15f %20.15f\n",
//                              id.c_str(), mglo, rmglo, Beta[mrki], Tdist[mrki], Se[mrki], Pval[mrki]);
//            assert(cx >= 0 && cx < LLEN);
//        }
//
//        // Collect numbers of markers to print in each task
//        int mp = M - n_rem;
//        int* mps = (int*) malloc(nranks * sizeof(int));
//        check_mpi(MPI_Allgather(&mp, 1, MPI_INTEGER, mps, 1, MPI_INTEGER, MPI_COMM_WORLD), __LINE__, __FILE__);
//
//        int ps = 0;
//        for (int i=0; i<rank; i++) { ps += mps[i]; }
//
//        MPI_Offset offset = size_t(ps) * size_t(LLEN-1);
//        check_mpi(MPI_File_write_at(*mlma_fh,
//                                    offset, todump, size_t(LLEN-1) * size_t(mp), MPI_CHAR, &status),
//                  __LINE__, __FILE__);
//
//        _mm_free(todump);
//
//        fflush(stdout);
//        double t_4 = MPI_Wtime();
//        MPI_Barrier(MPI_COMM_WORLD);
//        //if (rank >= 0)
//        //    printf("INFO   : intermediate time 4 = %.2f seconds.\n", t_4 - t_3);
//
//        MPI_Barrier(MPI_COMM_WORLD);
//
//        phen.close_prediction_files();
//
//        _mm_free(Beta);
//        _mm_free(Tdist);
//        _mm_free(Se);
//        _mm_free(Pval);
//        _mm_free(beta_sum);
//        _mm_free(beta_it);
//        _mm_free(g_k);
//        _mm_free(g);
//        _mm_free(y_k);
//    }
//
//    fflush(stdout);
//    double te = MPI_Wtime();
//    MPI_Barrier(MPI_COMM_WORLD);
//    if (rank == 0)
//        printf("INFO   : Time to compute the predictions: %.2f seconds.\n", te - ts);
//}

// EO: bim_file - I assume that the row number is the index
//
void Bayes::cross_bim_files() {

    if (rank == 0) {
        printf("INFO   : bim file:     %s\n", opt.get_bim_file().c_str());
        printf("INFO   : ref bim file: %s\n", opt.get_ref_bim_file().c_str());
    }
    std::ifstream in(opt.get_bim_file().c_str());
    if (!in) throw ("Error: can not open the file [" + opt.get_bim_file() + "] to read.");
    std::string   id, allele1, allele2;
    unsigned chr, physPos, idx = 0;
    float    genPos;
    while (in >> chr >> id >> genPos >> physPos >> allele1 >> allele2) {
        rsid.push_back(id);
    }
    in.close();
    int nrsid = rsid.size();
    if (rank == 0)
        printf("INFO   : found %d ids in bim file\n", nrsid);

    std::ifstream refin(opt.get_ref_bim_file().c_str());
    if (!refin) throw ("Error: can not open the file [" + opt.get_ref_bim_file() + "] to read.");
    idx = 0;
    while (refin >> chr >> id >> genPos >> physPos >> allele1 >> allele2) {
        m_refrsid[id] = idx++;
    }
    refin.close();
    if (rank == 0)
        printf("INFO   : found %d ids in reference bim file\n", idx);
}


// Multitrait
void Bayes::process() {

    //check_openmp();
    const int NPHEN = pmgr.get_phens().size(); // Number of phenotypes, pmgr object is created in BayesRR brr(opt, dims)

    /* 1.Initialization */
    pmgr.set_midx(); // Set marker ids
    pmgr.set_initial_sigmag(opt.get_ngroups(), NPHEN, Mm); // Set initial covariance matrix for all groups
    pmgr.set_pi_est(0.5, Mm); // Initialize vector of pi_priors for all markers
    pmgr.set_initial_sigmae(NPHEN, N); // Initialize covariance matrix for epsilon
    pmgr.set_initial_epsilon_mu(NPHEN); // setting initial mean as in intercept for phenotypes
    pmgr.set_acum(); // Setting vector to track number of times each marker has been processed
    pmgr.set_initial_betas(NPHEN); // Set initial betas as zeros for 1st iteration
    std::vector<std::vector<double>> psi = create_identity(NPHEN); // creation of initial identity matrix hyper


    //bool recv_update[nranks];

    // Create and initialize an output filestream object
    std::string rescov_fp = pmgr.get_outcsv_rescov_fp();
    std::ofstream file_rescov(rescov_fp);
    file_rescov << "it,00,01,10,11\n"; // Columns for output for covariance of residuals are iteration|00 elem of cov matrix|01|10|00

    std::string betascov_fp = pmgr.get_outcsv_betascov_fp();
    std::ofstream file_betascov(betascov_fp);
    string betas_cov_cols="it";
    for (int g = 0; g < G; ++g) {
        string gr = to_string(g);
        betas_cov_cols += ",00_" + gr + ",01_" + gr + ",10_" + gr + ",11_" + gr;
    };
    betas_cov_cols += "\n";
    file_betascov << betas_cov_cols; // Columns for output for covariance of betas are iteration and for each cell of matrix of each group

    std::string betas_fp = pmgr.get_outcsv_betas_fp();
    std::ofstream file_betas(betas_fp);
    file_betas << "it";
    for (int i = 0; i < 2000; ++i) {
        for (int j = 0; j < NPHEN; ++j) {
            file_betas << ","+to_string(i) + "_" + to_string(j);
        };
    };
    file_betas << "\n"; // Columns for output for betas are iteration and for marker for each pnehotype


    /* 2. Looping through iterations */
    for (int it = 1; it <= opt.get_iterations(); it++) {

        //double ts_it = MPI_Wtime();
        auto start = std::chrono::high_resolution_clock::now();
        int not_updated_num = 0;
        int updated_num = 0;
        

        /* 2.1. Add/subtract intercept from epsilon for each phenotype */
        int pidx = 0;
        for (auto& phen : pmgr.get_phens()) {
            phen.offset_epsilon(pmgr.get_epsilon_mu_phen(pidx));
            if (it == 1) {
                pmgr.set_sigmae_phen(pidx, phen.update_epsilon_sigma());
            };
            pidx++;
        };

        //if (it == 1) {
        //    pmgr.set_epsilon_mu(pmgr.get_epsilon_means());
        //}
        //else {
        //    pmgr.set_epsilon_mu(pmgr.sample_mvnorm_rng(pmgr.get_epsilon_means(), pmgr.multiply_matrix_scalar(pmgr.get_sigmae(),1/N), NPHEN));
        //};

        pmgr.set_epsilon_mu(pmgr.sample_mvnorm_rng(pmgr.get_epsilon_means(), pmgr.multiply_matrix_scalar(pmgr.get_sigmae(), 1.0/N), NPHEN));
        std::vector<double> t = pmgr.get_epsilon_mu();

        int pidx1 = 0;
        for (auto& phen : pmgr.get_phens()) {
            phen.offset_epsilon(-1 * pmgr.get_epsilon_mu_phen(pidx1));
            pidx1++;
        };

        // Inverse of covariance for epsilon for each iteration
        std::vector<std::vector<double>> sigmae_inv = pmgr.InvertMatrix(pmgr.get_sigmae(), NPHEN);


        if (opt.shuffle_markers()) {
            pmgr.shuffle_midx(opt.mimic_hydra());
        };

        /* 2.2. Looping through markers */
        for (int mrki=0; mrki<Mm; mrki++) {

            // To control the changes in beta of one marker, used for epsilon update
            std::vector<double> dbetas(NPHEN * 3); // dbeta, mave, msig for all phenotypes
            int mloc = 0;

            // If marker < M=task marker length
            if (mrki < M) {

                mloc = pmgr.get_marker_local_index(mrki); // Location index of marker
                const int mglo = S + mloc; // S=task marker start, mglo-global index of markers
                const int mgrp = get_marker_group(mglo); // Get group of the marker
                pmgr.add_marker_acum(mloc, 1); // tracking the number of times each marker has been processed

                // Get beta, sigma_group, sigma_epsilon and their inverses
                std::vector<double> beta = pmgr.get_marker_beta(mloc); // q-dim beta for concrete marker
                std::vector<std::vector<double>> sigmag = pmgr.get_sigmag_for_group(mgrp); // Covariance for group
                std::vector<std::vector<double>> sigmag_inv = pmgr.InvertMatrix(sigmag, NPHEN); // Inverse of covariance for group

                // Calculate denominator = covariance for sampling beta
                std::vector<std::vector<double>> denom = pmgr.sum_matrices(pmgr.multiply_matrix_scalar(sigmae_inv, (double)(N - 1.0)), sigmag_inv);
                std::vector<std::vector<double>> denom_inv = pmgr.InvertMatrix(denom, NPHEN);

                // Calculate numerator = Xj.T * (eps+XjBetaj) for the future mean. Xj.T * eps is calculated for each phenotype
                std::vector<double> num;
                for (auto& phen : pmgr.get_phens()) {
                    num.push_back(dot_product(mloc, phen.get_epsilon(), phen.get_marker_ave(mloc), phen.get_marker_sig(mloc)));
                };
                num = pmgr.sum_vectors(num, pmgr.multiply_vector_scalar(beta, (int)(N - 1)));
                num = pmgr.vec_matrix_product(num, sigmae_inv);
                num = pmgr.matrix_vec_product(denom_inv, num);


                // Calculate inclusion probability
                double f;
                f = 0.5 * pmgr.vec_vec_product(pmgr.vec_matrix_product(num, denom), num);
                if (f > 1000) {
                    f = 1000;
                };
                double pi0 = pmgr.get_pi_est(mloc);

                //std::cout << sqrt(pmgr.determinant(sigmag_inv, NPHEN)) << "\n";
                //std::cout << pow(pmgr.determinant(denom_inv, NPHEN), (double)(NPHEN / 2)) << "\n";
                //std::cout << exp(f) << "\n";

                double incl_prob = pi0 / (pi0 + (1 - pi0) * sqrt(pmgr.determinant(sigmag_inv, NPHEN)) * pow(pmgr.determinant(denom_inv, NPHEN), (double)(NPHEN / 2)) * exp(f));

                // Check if beta is sampled from MVN or not
                double prob = pmgr.sample_unif_rng();

                if (prob>incl_prob) {
                    //printf("In iteration %i for marker %i beta is updated. \n", it, mloc);
                    pmgr.set_marker_beta(mloc, pmgr.sample_mvnorm_rng(num, denom_inv, NPHEN));
                    pmgr.set_pi_for_marker(mloc, pmgr.sample_beta(1, 2));
                    updated_num++;
                }
                else {
                    //printf("In iteration %i for marker %i beta is not updated. \n", it, mloc);
                    pmgr.set_pi_for_marker(mloc, pmgr.sample_beta(2, 1));
                    not_updated_num++;
                };

                // Update epsilon if beta for the marker has been changed
                std::vector<double> beta_new = pmgr.get_marker_beta(mloc);
                int pheni = 0;
                for (auto& phen : pmgr.get_phens()) {
                    if (abs(beta[pheni] - beta_new[pheni]) > 0.0) {
                        dbetas[pheni * 3 + 0] = beta[pheni] - beta_new[pheni];
                        dbetas[pheni * 3 + 1] = phen.get_marker_ave(mloc);
                        dbetas[pheni * 3 + 2] = phen.get_marker_sig(mloc);
                    };
                    pheni++;
                };
                update_epsilon(mloc, dbetas); // Updating epsilon for each marker if beta !=0. Add Xj*(b_prev_j - b_new_j) to epsilon
            };

            //MPI_Barrier(MPI_COMM_WORLD);
            //double ts_sync = MPI_Wtime();

            // Collect information on which tasks need to share its just processed marker;
            // Tasks with M < Mm, share_mrk is false by default.
            //MPI_Allgather(&share_mrk,  1, MPI_C_BOOL,
            //              recv_update, 1, MPI_C_BOOL,
            //              MPI_COMM_WORLD);

            /*
            if (rank == 0) {
                printf("rank:share mrk %d?: ", mrki);
                for (int i=0; i<nranks; i++) {
                    printf("%d:%d ", i, recv_update[i] == true ? 1:0);
                }
                printf("\n");
            }
            fflush(stdout);
            */

            //int totbytes = 0;

            //int dis_bet[nranks], cnt_bet[nranks];
            //int dis_bed[nranks], cnt_bed[nranks];
            //int disp_bet = 0, disp_bed = 0;
            //for (int i=0; i<nranks; i++) {
            //    if (recv_update[i]) {
            //        cnt_bet[i] = NPHEN * 3;
            //        cnt_bed[i] = mbytes;
            //        totbytes += mbytes;
            //    } else {
            //        cnt_bet[i] = 0;
            //        cnt_bed[i] = 0;
            //    }
            //    dis_bet[i] = disp_bet;
            //    dis_bed[i] = disp_bed;
            //    //printf("mrki = %d, recv_update[%d] = %s %d:%d\n", mrki, i, recv_update[i] ? "T" : "F", dis_bet[i], cnt_bet[i]);
            //    disp_bet += cnt_bet[i];
            //    disp_bed += cnt_bed[i];
            //}

            //printf("rank %d bytes vs %d\n", totbytes, disp_bed);

            //double recv_dbetas[disp_bet];
            //MPI_Allgatherv(&dbetas, share_mrk ? NPHEN * 3 : 0, MPI_DOUBLE,
            //               recv_dbetas, cnt_bet, dis_bet, MPI_DOUBLE, MPI_COMM_WORLD);

            //unsigned char* recv_bed = (unsigned char*) _mm_malloc(disp_bed, 32);
            //check_malloc(recv_bed, __LINE__, __FILE__);
            //MPI_Allgatherv(&bed_data[mloc * mbytes], share_mrk ? mbytes : 0, MPI_UNSIGNED_CHAR,
            //               recv_bed, cnt_bed, dis_bed, MPI_UNSIGNED_CHAR, MPI_COMM_WORLD);

            //MPI_Barrier(MPI_COMM_WORLD);
            //double te_sync = MPI_Wtime();
            //t_it_sync += te_sync - ts_sync;

            //_mm_free(recv_bed);
        }; // End marker loop

        //int pheni = -1;
        //for (auto& phen : pmgr.get_phens()) {
        //    pheni++;

        //    phen.reset_beta_sqn_to_zero();
        //    for (int i = 0; i < M; i++) {
        //        phen.increment_beta_sqn(get_marker_group(S + i), phen.get_marker_beta(i) * phen.get_marker_beta(i)); // Adding value to position in sequence
        //    }

        //    double* beta_sqn_sum = (double*)malloc(G * sizeof(double));
        //    check_mpi(MPI_Allreduce(phen.get_beta_sqn()->data(),
        //        beta_sqn_sum,
        //        G,
        //        MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD), __LINE__, __FILE__);
        //    phen.set_beta_sqn(beta_sqn_sum);
        //    free(beta_sqn_sum);

        //    int* cass_sum = (int*)malloc(G * K * sizeof(int));
        //    check_mpi(MPI_Allreduce(phen.get_cass(),
        //        cass_sum,
        //        G * K,
        //        MPI_INT,
        //        MPI_SUM,
        //        MPI_COMM_WORLD), __LINE__, __FILE__);
        //    phen.set_cass(cass_sum);
        //    free(cass_sum);
        //}

        /* 2.3. Update covariance for each group */
        for (int i=0; i<G; i++) {

            // Skip empty groups - WHAT IS IT
            //if (mtotgrp.at(i) == 0)
            //    continue;

            //pmgr.set_m0_for_group(i, mtotgrp.at(i) - pmgr.get_cass_for_group(i, 0));
            // Skip groups with m0 being null or empty cass (adaV in action)  - WHAT IS IT
            //if (phen.get_m0_for_group(i) == 0 || phen.get_cass_sum_for_group(i) == 0) {
            //    phen.set_sigmag_for_group(i, 0.0);
            //    continue;
            //}
            std::vector<std::vector<double>> betas_group;
            for (int j = 0; j < Mm; j++) {
                if (group_index[j] == i) {
                    std::vector<double> b = pmgr.get_marker_beta(j);
                    betas_group.push_back(b);
                };
            };

            // Update covariance for each group from InvWishart
            pmgr.set_sigmag_for_group(i, pmgr.sample_inv_wishart_rng(NPHEN, mtotgrp[i]+NPHEN, pmgr.sum_matrices(psi, pmgr.mul_matrices(betas_group, betas_group))));
        }

        // Broadcast sigmaG of rank 0
        //check_mpi(MPI_Bcast(phen.get_sigmag()->data(), phen.get_sigmag()->size(), MPI_DOUBLE, 0, MPI_COMM_WORLD), __LINE__, __FILE__);
        //double e_sqn = phen.epsilon_sumsqr();

        /* 2.4. Update covariance for epsilon */
        pmgr.set_sigmae(pmgr.sample_inv_wishart_rng(NPHEN, N + NPHEN, pmgr.sum_matrices(psi, eps_eps_product(NPHEN))));

        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
        cout << "The duration of iteration "<<it<<" was " << duration.count() << endl;
        printf("For iteration %i number of upated markers %i and of not updated markers %i.\n",it,updated_num,not_updated_num);

        /* 2.5. Push data to output files */
        if (it >= 1000) {
            // Write row to residuals covariance file
            file_rescov << it;
            std::vector<std::vector<double>> sigmae_new = pmgr.get_sigmae();
            for (int j = 0; j < NPHEN; j++)
            {
                for (int k = 0; k < NPHEN; k++)
                {
                    file_rescov << ","; 
                    file_rescov << sigmae_new[j][k];
                };
            };
            file_rescov << "\n";

            // Write row to betas covariance file
            file_betascov << it;
            for (int gr = 0; gr < G; ++gr) {
                std::vector<std::vector<double>> sigmag_gr_new = pmgr.get_sigmag_for_group(gr);
                for (int j = 0; j < NPHEN; j++)
                {
                    for (int k = 0; k < NPHEN; k++)
                    {
                        file_betascov << ",";
                        file_betascov << sigmag_gr_new[j][k];
                    };
                };
            };
            file_betascov << "\n";

            // Write row to betas file
            file_betas << it;
            std::vector<std::vector<double>> betas = pmgr.get_betas();
            for (int i = 0; i < 2000; ++i) {
                for (int j = 0; j < NPHEN; ++j) {
                    file_betas << ",";
                    file_betas << betas[i][j];
                };
            };
            file_betas << "\n";

        };


        //double sigmae_r0 = pmgr.get_sigmae();
        //check_mpi(MPI_Bcast(&sigmae_r0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD), __LINE__, __FILE__);
        //phen.set_pmgr(sigmae_r0);
        //if (rank % 10 == 0) {
        //    printf("RESULT : i:%d r:%d p:%d  sum sigmaG = %20.15f  sigmaE = %20.15f\n", it, rank, pheni, phen.get_sigmag_sum(), phen.get_sigmae());
        //}

        // Broadcast pi_est from rank 0
        //for (int i=0; i<G; i++) {
        //    check_mpi(MPI_Bcast(phen.get_pi_est()->at(i).data(), phen.get_pi_est()->at(i).size(), MPI_DOUBLE, 0, MPI_COMM_WORLD), __LINE__, __FILE__);
        //}

        //double te_it = MPI_Wtime();
        //if (rank == 0)
        //    printf("RESULT : It %d  total proc time = %7.3f sec, with sync time = %7.3f\n", it, te_it - ts_it, t_it_sync);


        // Write output files
        //if (it % opt.get_output_thin_rate() == 0) {
        //    const unsigned nthinned = it / opt.get_output_thin_rate() - 1;
        //    for (auto& phen : pmgr.get_phens()) {
        //        if (rank == 0) {
        //            //phen.print_pi_est();
        //            write_ofile_csv(*(phen.get_outcsv_fh()), it,  phen.get_sigmag(), phen.get_sigmae(), phen.get_m0_sum(), nthinned, phen.get_pi_est());
        //        }
        //        write_ofile_h1(*(phen.get_outbet_fh()), rank, Mt, it, nthinned, S, M, phen.get_betas().data(), MPI_DOUBLE);
        //        write_ofile_h1(*(phen.get_outcpn_fh()), rank, Mt, it, nthinned, S, M, phen.get_comp().data(),  MPI_INTEGER);
        //    }
        //}

    } // End iteration loop

    file_rescov.close();
    file_betascov.close();
    file_betas.close();

    //MPI_Barrier(MPI_COMM_WORLD);

    //for (auto& phen : pmgr.get_phens())
    //    phen.close_output_files();

} // End process function

std::vector<std::vector<double>> Bayes::create_identity(const int q) {
    std::vector<std::vector<double>> iden(q, std::vector<double>(q));
    for (int j = 0; j < q; j++) {
        for (int k = 0; k < q; k++) {
            if (k == j) {
                iden[k][j] = 1.0;
            };
        };
    };

    return iden;
}


void Bayes::update_epsilon(const int mloc, const std::vector<double> dbetas) {

    int pheni = 0;
    for (auto& phen : pmgr.get_phens()) {
        if (dbetas[pheni*3] != 0.0) {
            phen.update_epsilon(&dbetas[pheni*3], &bed_data[mloc * mbytes]);
        };
        pheni += 1;
    };
}


std::vector<std::vector<double>> Bayes::eps_eps_product(const int q) {

    std::vector<std::vector<double>> eps_prod(q, std::vector<double>(q));
    unsigned char* bed = &bed_data[0];

    for (int q1 = 0; q1 < q; q1++) {
        for (int q2 = 0; q2 <= q1; q2++) {

            int pheni1 = 0;
            for (auto& phen1 : pmgr.get_phens()) {
                if (pheni1 == q1) {
                    double* eps_phen1 = phen1.get_epsilon();
                    int pheni2 = 0;
                    for (auto& phen2 : pmgr.get_phens()) {
                        if (pheni2 == q2) {
                            double* eps_phen2 = phen2.get_epsilon();

                            for (int i = 0; i < mbytes; i++) {
                                for (int j = 0; j < 4; j++) {
                                    //std::cout << q1 << "\n" << q2 << "\n" << i << "\n" << j << "\n";
                                    //std::cout << eps_phen1[i * 4 + j]<<"\n"<< eps_phen2[i * 4 + j] << "\n";

                                    eps_prod[q1][q2] += eps_phen1[i * 4 + j] * eps_phen2[i * 4 + j];
                                    if (q1 != q2) {
                                        eps_prod[q2][q1] += eps_phen1[i * 4 + j] * eps_phen2[i * 4 + j];
                                    };
                                };
                            };

                        };
                        pheni2++;
                    };
                };
                pheni1++;
            };


        };
    };

    return eps_prod;
}


//std::vector<std::vector<double>> Bayes::eps_eps_product1(const int q) {
//
//    std::vector<std::vector<double>> eps_prod(q, std::vector<double>(q));
//    unsigned char* bed = &bed_data[0];
//
//    //int pheni = 0;
//    //for (auto& phen : pmgr.get_phens()) {
//    //    double* eps_phen = phen.get_epsilon();
//    //    for (int i = 0; i < mbytes; i++) {
//    //        for (int j = 0; j < 4; j++) {
//    //            eps_prod[pheni][pheni] += eps_phen[i * 4 + j] * eps_phen[i * 4 + j];
//    //        };
//    //    };
//    //    pheni++;
//    //};
//
//    int pheni1 = 0;
//    for (auto& phen1 : pmgr.get_phens()) {
//        double* eps_phen1 = phen1.get_epsilon();
//        int pheni2 = 0;
//        for (auto& phen2 : pmgr.get_phens()) {
//            double* eps_phen2 = phen2.get_epsilon();
//
//
//            for (int i = 0; i < mbytes; i++) {
//                for (int j = 0; j < 4; j++) {
//                    eps_prod[pheni1][pheni2] += eps_phen1[i * 4 + j] * eps_phen2[i * 4 + j];
//                };
//            };
//            pheni2++;
//        };
//        pheni1++;
//
//    };
//
//    return eps_prod;
//}



//void Bayes::update_epsilon(const int* counts, const double* dbetas, const unsigned char* bed) {
//
//    const int NPHEN =  pmgr.get_phens().size();
//
//    int cnt_tot = 0;
//    int bedi = 0;
//    for (int i=0; i<nranks; i++) {
//        int cnt = counts[i];
//        if (cnt == 0) continue;
//        assert(cnt == NPHEN * 3);
//        int pheni = 0;
//        for (auto& phen : pmgr.get_phens()) {
//            //printf("r:%d p:%d - dbeta = %20.15f\n", rank, pheni, dbetas[cnt_tot + pheni * 3]);
//            if (dbetas[cnt_tot + pheni * 3] != 0.0)
//                phen.update_epsilon(&dbetas[cnt_tot + pheni * 3], &bed[bedi * mbytes]);
//            //phen.update_epsilon_sum();
//            //printf(" - r:%d p:%d   epssum = %.17g, cnt_tot = %d\n", i, pheni, phen.get_epsilon_sum(), cnt_tot);
//            pheni += 1;
//        }
//        cnt_tot += cnt;
//        bedi    += 1;
//    }
//    fflush(stdout);
//}


// Product for one phenotype epsilon
double Bayes::dot_product(const int mloc, double* eps_phen, const double mu_marker, const double sigma_marker) {

    unsigned char* bed = &bed_data[mloc * mbytes];

//#ifdef MANVEC
//    __m256d luta, lutb; //, lutna;
//    __m512d lutab, p42;
//    __m256d p4 = _mm256_set1_pd(0.0);
//    __m256d suma = _mm256_set1_pd(0.0);
//    __m256d sumb = _mm256_set1_pd(0.0);
//    __m512d sum42 = _mm512_set1_pd(0.0);
//
//#ifdef _OPENMP
//    //#pragma omp parallel for schedule(static) reduction(addpd4:suma,sumb)
//#pragma omp parallel for schedule(static) reduction(addpd8:sum42)
//#endif
//    for (int j = 0; j < mbytes; j++) {
//        //luta  = _mm256_load_pd(&dotp_lut_a[bed[j] * 4]);
//        //lutb  = _mm256_load_pd(&dotp_lut_b[bed[j] * 4]);
//        //luta  = _mm256_load_pd(&dotp_lut_ab[bed[j] * 8]);
//        //lutb  = _mm256_load_pd(&dotp_lut_ab[bed[j] * 8 + 4]);
//        lutab = _mm512_load_pd(&dotp_lut_ab[bed[j] * 8]);
//        //p4    = _mm256_load_pd(&phen[j * 4]);
//        p42 = _mm512_broadcast_f64x4(_mm256_load_pd(&phen[j * 4]));
//        ////lutna = _mm256_load_pd(&na_lut[mask4[j] * 4]); // phen = 0.0 on NAs!
//        //luta  = _mm256_mul_pd(luta, p4);
//        //lutb  = _mm256_mul_pd(lutb, p4);
//        p42 = _mm512_mul_pd(p42, lutab);
//        //suma  = _mm256_add_pd(suma, luta);
//        //sumb  = _mm256_add_pd(sumb, lutb);
//        sum42 = _mm512_add_pd(sum42, p42);
//    }
//
//    //return sigma_inv *
//    //    (suma[0] + suma[1] + suma[2] + suma[3] - mu * (sumb[0] + sumb[1] + sumb[2] + sumb[3]));
//    return sigma_inv *
//        (sum42[0] + sum42[1] + sum42[2] + sum42[3] - mu * (sum42[4] + sum42[5] + sum42[6] + sum42[7]));
//
//#else

    double dpa = 0.0;
    double dpb = 0.0;

//#ifdef _OPENMP
//    //#pragma omp parallel for schedule(static) reduction(addpd4:suma,sumb)
//#pragma omp parallel for schedule(static) reduction(+:dpa,dpb)
//#endif
    for (int i = 0; i < mbytes; i++) {
//#ifdef _OPENMP
//#pragma omp simd aligned(dotp_lut_a,dotp_lut_b,phen:32)
//#endif
        for (int j = 0; j < 4; j++) {
            dpa += dotp_lut_a[bed[i] * 4 + j] * eps_phen[i * 4 + j];
            dpb += dotp_lut_b[bed[i] * 4 + j] * eps_phen[i * 4 + j];
            //std::cout << "dotp_lut_a" << "\n"<<dotp_lut_a[bed[i] * 4 + j] << "\n";
            //std::cout << "dotp_lut_b" << "\n"<< dotp_lut_b[bed[i] * 4 + j] << "\n";
            //std::cout << "eps_phen" << "\n"<<eps_phen[i * 4 + j] << "\n";
        }
    }

    return (dpa - mu_marker * dpb) * sigma_marker; // sigma is already inversed

//#endif

}


// Setup processing: load input files and define MPI task workload
void Bayes::setup_processing() {

    mbytes = (N %  4) ? (size_t) N /  4 + 1 : (size_t) N /  4;
    load_genotype();

    pmgr.read_phen_files(opt, get_N(), get_M());
    pmgr.set_output_filenames(opt.get_out_dir()); // Setting output files for covariance of epsilon, covariances for betas.

    check_processing_setup();

    //if (rank == 0)
    //    printf("INFO   : output directory: %s\n", opt.get_out_dir().c_str());

    //MPI_Barrier(MPI_COMM_WORLD);
    //double ts = MPI_Wtime();
    pmgr.compute_markers_statistics(bed_data, get_N(), get_M(), mbytes);
    //pmgr.display_markers_statistics(get_M());
    //MPI_Barrier(MPI_COMM_WORLD);
    //double te = MPI_Wtime();
    //if (rank == 0)
    //    printf("INFO   : Time to compute the markers' statistics: %.2f seconds.\n", te - ts);

    //if (opt.predict()) return;

    //for (auto& phen : pmgr.get_phens()) {
    //    phen.set_prng_m((unsigned int)(opt.get_seed() + (rank + 0)));
    //    if (opt.mimic_hydra()) {
    //        phen.set_prng_d((unsigned int)(opt.get_seed() + (rank + 0) * 1000));
    //    } else {
    //        phen.set_prng_d((unsigned int)(opt.get_seed() + (rank + 1) * 1000));
    //    }
    //}

    pmgr.set_prng_m((unsigned int)(opt.get_seed() + (rank + 0)));
    if (opt.mimic_hydra()) {
        pmgr.set_prng_d((unsigned int)(opt.get_seed() + (rank + 0) * 1000));
    }
    else {
        pmgr.set_prng_d((unsigned int)(opt.get_seed() + (rank + 1) * 1000));
    }

    read_group_index_file(opt.get_group_index_file());

    for (int i = 0; i < Mt; i++) {
        mtotgrp.at(get_marker_group(i)) += 1;
    };

    //for (int i=0; i<G; i++)
    //    printf("mtotgrp at %d = %d\n", i, mtotgrp.at(i));
}


//void Bayes::check_openmp() {
//#ifdef _OPENMP
//#pragma omp parallel
//    {
//        if (rank == 0) {
//            int nt = omp_get_num_threads();
//            if (omp_get_thread_num() == 0)
//                printf("INFO   : OMP parallel regions will use %d thread(s)\n", nt);
//        }
//    }
//#else
//    printf("WARNING: no OpenMP support!\n");
//#endif
//}

void Bayes::read_group_index_file(const std::string& file) {

    std::ifstream infile(file.c_str());
    if (! infile)
        throw ("Error: can not open the group file [" + file + "] to read. Use the --group-index-file option!");

    //if (rank == 0)
    //    std::cout << "INFO   : Reading groups from " + file + "." << std::endl;

    std::string label;
    int group;

    group_index.clear();

    while (infile >> label >> group) {
        //std::cout << label << " - " << group << std::endl;
        if (group > G) {
            printf("FATAL  : group index file contains a value that exceeds the number of groups given in group mixture file.\n");
            printf("       : check the consistency between your group index and mixture input files.\n");
            exit(1);
        }
        group_index.push_back(group);
    }
}


void Bayes::check_processing_setup() {
    for (const auto& phen : pmgr.get_phens()) {
        const int Np = phen.get_nas() + phen.get_nonas();
        if (Np != N) {
            std::cout << "Fatal: N = " << N << " while phen file " << phen.get_filepath() << " has " << Np << " individuals!" << std::endl;
            exit(1);
        }
    }
}


void Bayes::load_genotype() {

    //double ts = MPI_Wtime();

    //MPI_File bedfh;
    const std::string bedfp = opt.get_bed_file();
    //check_mpi(MPI_File_open(MPI_COMM_WORLD, bedfp.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &bedfh),  __LINE__, __FILE__);

    const size_t size_bytes = size_t(M) * size_t(mbytes) * sizeof(unsigned char);

    bed_data = (unsigned char*)_mm_malloc(size_bytes, 64);

    //ifstream bedfs(bedfp, ios::binary);
    ifstream bedfs;
    bedfs.open(bedfp.c_str(), ios::in | ios::binary);
    bedfs.seekg(3, ios::beg);
    bedfs.read((char*)bed_data, size_bytes);
    bedfs.close();


    //check_malloc(bed_data, __LINE__, __FILE__);
    //printf("INFO   : rank %4d has allocated %zu bytes (%.3f GB) for raw data.\n", rank, size_bytes, double(size_bytes) / 1.0E9);

    // Offset to section of bed file to be processed by task
    //MPI_Offset offset = size_t(3) + size_t(S) * size_t(mbytes) * sizeof(unsigned char);

    // Gather the sizes to determine common number of reads
    size_t max_size_bytes = 0;
    //check_mpi(MPI_Allreduce(&size_bytes, &max_size_bytes, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD), __LINE__, __FILE__);

    //const int NREADS = check_int_overflow(size_t(ceil(double(max_size_bytes)/double(INT_MAX/2))), __LINE__, __FILE__);
    size_t bytes = 0;
    //mpi_file_read_at_all <unsigned char*> (size_bytes, offset, bedfh, MPI_UNSIGNED_CHAR, NREADS, bed_data, bytes);
    //@@@MPI_Barrier(MPI_COMM_WORLD);

    //check_mpi(MPI_File_close(&bedfh), __LINE__, __FILE__);

    //double te = MPI_Wtime();
    //MPI_Barrier(MPI_COMM_WORLD);
    //if (rank >= 0)
    //    printf("INFO   : time to load genotype data = %.2f seconds.\n", te - ts);

}


//void Bayes::set_block_of_markers() {
//
//    const int modu = Mt % nranks;
//    const int size = Mt / nranks;
//
//    Mm = Mt % nranks != 0 ? size + 1 : size;
//
//    int len[nranks], start[nranks];
//    int cum = 0;
//    for (int i=0; i<nranks; i++) {
//        len[i]  = i < modu ? size + 1 : size;
//        start[i] = cum;
//        cum += len[i];
//    }
//    //printf("cum %d vs %d  Mm = %d\n", cum, Mt, Mm);
//    assert(cum == Mt);
//
//    M = len[rank];
//    S = start[rank];
//
//    printf("INFO   : rank %4d has %d markers over tot Mt = %d, max Mm = %d, starting at S = %d\n", rank, M, Mt, Mm, S);
//    //@todo: mpi check sum over tasks == Mt
//}

void Bayes::set_block_of_markers() {

    //const int modu = Mt % nranks;
    //const int size = Mt / nranks;

    //Mm = Mt % nranks != 0 ? size + 1 : size;

    //int len[nranks], start[nranks];
    //int cum = 0;
    //for (int i=0; i<nranks; i++) {
    //    len[i]  = i < modu ? size + 1 : size;
    //    start[i] = cum;
    //    cum += len[i];
    //}
    ////printf("cum %d vs %d  Mm = %d\n", cum, Mt, Mm);
    //assert(cum == Mt);

    //M = len[rank];
    //S = start[rank];

    //printf("INFO   : rank %4d has %d markers over tot Mt = %d, max Mm = %d, starting at S = %d\n", rank, M, Mt, Mm, S);
    //@todo: mpi check sum over tasks == Mt
    M = Mt;
    Mm = Mt;
}

//void Bayes::print_cva() {
//    printf("INFO   : mixtures for all groups:\n");
//    for (int i=0; i<G; i++) {
//        printf("         grp %2d: ", i);
//        for (int j=0; j<K; j++) {
//            printf("%7.5f ", cva[i][j]);
//        }
//        printf("\n");
//    }
//}
//
//void Bayes::print_cvai() {
//    printf("INFO   : inverse mixtures for all groups:\n");
//    for (int i=0; i<G; i++) {
//        printf("         grp %2d: ", i);
//        for (int j=0; j<K; j++) {
//            printf("%10.3f ", cvai[i][j]);
//        }
//        printf("\n");
//    }
//}
