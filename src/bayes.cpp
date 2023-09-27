#include <iostream>
#include <fstream>
#include <iterator>
#include <limits.h>
#include <cmath>
#include <immintrin.h>
#include <omp.h>
#include "bayes.hpp"
#include "utilities.hpp"
#include "dotp_lut.hpp"
#include "na_lut.hpp"
#include "xfiles.hpp"
#include <boost/math/special_functions/gamma.hpp>
#include <regex>

void Bayes::predict() {

    MPI_Barrier(MPI_COMM_WORLD);
    double ts = MPI_Wtime();

    check_openmp();

    MPI_Status status;
    MPI_Offset file_size = 0;
    MPI_File*  fh;

    //cross_bim_files();

    int C = opt.get_cov_num();

    int pidx = 0;
    for (auto& phen : pmgr.get_phens()) {
        pidx += 1;

        phen.delete_output_prediction_files();
        phen.open_prediction_files();

        const std::vector<unsigned char> mask4 = phen.get_mask4();
        const int im4 = phen.get_im4();

        fh = phen.get_inbet_fh();
        check_mpi(MPI_File_get_size(*fh, &file_size), __LINE__, __FILE__);
        printf("file_size = %u B\n", file_size);
        std::cout << "file name: " << phen.get_inbet_fp() << std::endl;

        // First element of the .bet is the total number of processed markers
        // Then: iteration (uint) beta (double) for all markers
        uint Mtot_ = 0;
        MPI_Offset betoff = size_t(0);
        check_mpi(MPI_File_read_at_all(*fh, betoff, &Mtot_, 1, MPI_UNSIGNED, &status), __LINE__, __FILE__);
        if (Mtot_ != Mt) {
            printf("Mismatch between expected and Mtot read from .bet file: %d vs %d\n", Mt, Mtot_);
            MPI_Abort(MPI_COMM_WORLD, 1);
        }

        assert((file_size - sizeof(uint)) % (Mtot_ * sizeof(double) + sizeof(uint)) == 0);
        uint niter = (file_size - sizeof(uint)) / (Mtot_ * sizeof(double) + sizeof(uint));
        if (rank == 0)
            printf("INFO   : Number of recorded iterations in .bet file %d: %u\n", pidx-1, niter);

        double* beta_sum = (double*) _mm_malloc(size_t(Mtot_) * sizeof(double), 32);
        check_malloc(beta_sum, __LINE__, __FILE__);
        for (int i=0; i<Mtot_; i++) beta_sum[i] = 0.0;

        double* beta_it = (double*) _mm_malloc(size_t(Mtot_) * sizeof(double), 32);
        check_malloc(beta_it, __LINE__, __FILE__);

        uint start_iter = 100;
        //EO: use this one to speed up testing (avoids to read the entire bet history)
        //if (niter > 3) start_iter = niter - 3;

        for (uint i=start_iter; i<niter; i++) {
            betoff
                = sizeof(uint) // Mtot
                + (sizeof(uint) + size_t(Mtot_) * sizeof(double)) * size_t(i)
                + sizeof(uint);
            check_mpi(MPI_File_read_at_all(*fh, betoff, beta_it, Mtot_, MPI_DOUBLE, &status), __LINE__, __FILE__);
            for (int j=0; j<Mtot_;j++)
                beta_sum[j] += beta_it[j];
        }

        for (int j=0; j<Mtot_;j++)
            beta_sum[j] /= double(niter);
        /*
        if(C > 0){
            phen.load_cov_deltas();
            phen.set_Z(Z);
            phen.avg_deltas_it(niter);
        }
        std::vector<double>* deltas = phen.get_deltas();

        double* c = (double*) _mm_malloc(size_t(N) * sizeof(double), 32);
        check_malloc(c, __LINE__, __FILE__);
        
        for (int i=0; i<N; i++){
            c[i] = 0.0;
            for(int covi = 0; covi<C; covi++){
                c[i] += deltas->at(covi) * Z[i][covi];
            }
        }*/
        
        fflush(stdout);
        double t_1 = MPI_Wtime();
        MPI_Barrier(MPI_COMM_WORLD);
        //if (rank >= 0)
        //    printf("INFO   : intermediate time 1 = %.2f seconds.\n", t_1 - ts);


        double* g_k = (double*) _mm_malloc(size_t(N) * sizeof(double), 32);
        check_malloc(g_k, __LINE__, __FILE__);
        for (int i=0; i<N; i++) g_k[i] = 0.0;

#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
        for (int mrki=0; mrki<M; mrki++) {

            // Get rsid from current
            int mglo = S + mrki;
            //std::string id = rsid.at(mglo);

            // Skip markers with no corresponding rsid in reference bim file
            //if (m_refrsid.find(id) == m_refrsid.end()) {
                //printf("%d -> %d = %s not found in reference bim\n", mrki, mglo, id.c_str());
            //    continue;
            //}

            //int rmglo = m_refrsid.find(id)->second; // global marker index in reference bed

            size_t methix = size_t(mrki) * size_t(N);
            const double* methm = &meth_data[methix];

            double mave = phen.get_marker_ave(mrki);
            double msig = phen.get_marker_sig(mrki);

            for (int i = 0; i < N; i++) {

                    double val = (methm[i] - mave);
#ifdef _OPENMP
#pragma omp atomic update
#endif
                    g_k[i] += val * beta_sum[mglo];

            }
        }

        fflush(stdout);
        double t_2 = MPI_Wtime();
        MPI_Barrier(MPI_COMM_WORLD);
        //if (rank >= 0)
        //    printf("INFO   : intermediate time 2 = %.2f seconds.\n", t_2 - t_1);

        double* g = (double*) _mm_malloc(size_t(N) * sizeof(double), 32);
        check_malloc(g, __LINE__, __FILE__);

        check_mpi(MPI_Allreduce(g_k, g, N, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD), __LINE__, __FILE__);

        double* y_k = (double*) _mm_malloc(size_t(N) * sizeof(double), 32);
        check_malloc(y_k, __LINE__, __FILE__);

        phen.get_centered_and_scaled_y(y_k);

        std::ofstream yest_stream;
        yest_stream.open(phen.get_outyest_fp());
        for (int i=0; i<N; i++){
            double z = g[i];// + c[i];
            yest_stream << z << std::endl;
        }
        yest_stream.close();

        //EO: already done when computing original epsilon when reading .phen
	//phen.set_nas_to_zero(y_k, im4*4);
        // NAs are 0.0 in all, so should be safe
        for (int i=0; i<N; i++)
            y_k[i] -= (g[i] - g_k[i]);

        double sigma = 0.0;
        for (int i=0; i<N; i++)
            sigma += y_k[i] * y_k[i];
        sigma /= phen.get_nonas();
        //printf("### r: %d sigma = %20.15f\n", rank, sigma);
        /*
        MPI_File* mlma_fh = phen.get_outmlma_fh();

        double* Beta  = (double*) _mm_malloc(size_t(M) * sizeof(double), 32);
        check_malloc(Beta, __LINE__, __FILE__);
        double* Tdist = (double*) _mm_malloc(size_t(M) * sizeof(double), 32);
        check_malloc(Tdist, __LINE__, __FILE__);
        double* Se    = (double*) _mm_malloc(size_t(M) * sizeof(double), 32);
        check_malloc(Se, __LINE__, __FILE__);
        double* Pval  = (double*) _mm_malloc(size_t(M) * sizeof(double), 32);
        check_malloc(Pval, __LINE__, __FILE__);


#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (int mrki=0; mrki<M; mrki++) {

            // Global index in current .bim
            int mglo = S + mrki;
            std::string id = rsid.at(mglo);

            // Skip markers with no corresponding rsid in reference bim file
            if (m_refrsid.find(id) == m_refrsid.end()) {
                //printf("%d -> %d = %s not found in reference bim\n", mrki, mglo, id.c_str());
                continue;
            }

            // Global marker index in reference bed
            int rmglo = m_refrsid.find(id)->second;

            size_t methix = size_t(mrki) * size_t(N);
            const double* methm = &meth_data[methix];

            double mave = phen.get_marker_ave(mrki);
            double msig = phen.get_marker_sig(mrki);

            double xtx = 0.0;
            double xty = 0.0;
            double chk = 0.0;
            for (int i=0; i < N; i++) {
                
                    double val = methm[i] - mave;
                    xtx += val * val;
                    xty += val * y_k[i];
                
            }

            double beta  = xty / xtx;
            double tdist = xty / sqrt(sigma * xtx);
            double se    = beta / tdist;
            double pval  = 1.0 - boost::math::gamma_p(0.5, tdist * tdist * 0.5);

            Beta[mrki]  = beta;
            Tdist[mrki] = tdist;
            Se[mrki]    = se;
            Pval[mrki]  = pval;

            //if (mrki < 5) {
            //    printf("r: %3d  m: %8d %8d (%5s) %8d  beta=%20.15f  se=%20.15f  tdist=%20.15f pval=%20.15f   xtx=%10.1f xty=%15.6f\n", rank, mrki, mglo, id.c_str(), rmglo, beta, se, tdist, pval, xtx, xty);
            //}
        }

        fflush(stdout);
        double t_3 = MPI_Wtime();
        MPI_Barrier(MPI_COMM_WORLD);
        //if (rank >= 0)
        //    printf("INFO   : intermediate time 3 = %.2f seconds.\n", t_3 - t_2);
    
        const int LLEN = 123 + 1;
        char* todump  = (char*) _mm_malloc(size_t(LLEN) * size_t(M) * sizeof(char), 32);
        check_malloc(todump, __LINE__, __FILE__);

        int n_rem = 0;
        for (int mrki=0; mrki<M; mrki++) {
            int mglo = S + mrki;
            std::string id = rsid.at(mglo);
            if (m_refrsid.find(id) == m_refrsid.end()) {
                //|| id.compare("rs12562034") == 0m|| id.compare("rs188466450") == 0)
                printf("WARNING: marker id %s excluded -- no match\n", id.c_str());
                n_rem++;
                continue;
            }
            int rmglo = m_refrsid.find(id)->second;
            int cx = snprintf(&todump[(mrki - n_rem) * (LLEN - 1)], LLEN, "%20s %8d %8d %20.15f %20.15f %20.15f %20.15f\n",
                              id.c_str(), mglo, rmglo, Beta[mrki], Tdist[mrki], Se[mrki], Pval[mrki]);
            assert(cx >= 0 && cx < LLEN);
        }

        // Collect numbers of markers to print in each task
        int mp = M - n_rem;
        int* mps = (int*) malloc(nranks * sizeof(int));
        check_mpi(MPI_Allgather(&mp, 1, MPI_INTEGER, mps, 1, MPI_INTEGER, MPI_COMM_WORLD), __LINE__, __FILE__);

        int ps = 0;
        for (int i=0; i<rank; i++) { ps += mps[i]; }

        MPI_Offset offset = size_t(ps) * size_t(LLEN-1);
        check_mpi(MPI_File_write_at(*mlma_fh,
                                    offset, todump, size_t(LLEN-1) * size_t(mp), MPI_CHAR, &status),
                  __LINE__, __FILE__);

        _mm_free(todump);

        fflush(stdout);
        double t_4 = MPI_Wtime();
        MPI_Barrier(MPI_COMM_WORLD);
        //if (rank >= 0)
        //    printf("INFO   : intermediate time 4 = %.2f seconds.\n", t_4 - t_3);

        MPI_Barrier(MPI_COMM_WORLD);

        phen.close_prediction_files();

        _mm_free(Beta);
        _mm_free(Tdist);
        _mm_free(Se);
        _mm_free(Pval);
        _mm_free(beta_sum);
        _mm_free(beta_it);
        _mm_free(g_k);
        _mm_free(g);
        _mm_free(y_k);
        */
    }
    
    fflush(stdout);
    double te = MPI_Wtime();
    MPI_Barrier(MPI_COMM_WORLD);
    if (rank == 0)
        printf("INFO   : Time to compute the predictions: %.2f seconds.\n", te - ts);
}

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

void Bayes::process() {

    double ts_overall = MPI_Wtime();

    check_openmp();

    for (auto& phen : pmgr.get_phens()) {
        phen.delete_output_files();
        phen.open_output_files();
        phen.set_midx();
        for (int i=0; i<opt.get_ngroups(); i++) {
            phen.set_sigmag_for_group(i, phen.sample_beta_rng(1.0, 1.0));
            if (mtotgrp.at(i) == 0)
                phen.set_sigmag_for_group(i, 0.0);
            //printf("sample sigmag[%d] = %20.15f\n", i, phen.get_sigmag_for_group(i));
        }
        check_mpi(MPI_Bcast(phen.get_sigmag()->data(), phen.get_sigmag()->size(), MPI_DOUBLE, 0, MPI_COMM_WORLD), __LINE__, __FILE__);
        phen.set_pi_est(pi_prior);
        //printf("sample sigmag[0] = %20.15f\n", phen.get_sigmag()->at(0));
    }

    // number of covariates
    int C = opt.get_cov_num();

    const int NPHEN =  pmgr.get_phens().size();
    bool recv_update[nranks];

    for (unsigned int it = 1; it <= opt.get_iterations(); it++) {

        double ts_it = MPI_Wtime();

        if (rank == 0)
            printf("\n\n@@@ ITERATION %5d\n", it);

        int pidx = 0;
        for (auto& phen : pmgr.get_phens()) {
            pidx += 1;
            //printf("phen %d has mu = %20.15f\n", pidx, phen.get_mu());
            phen.offset_epsilon(phen.get_mu());
            //phen.update_epsilon_sum();
            if (it == 1) {
                phen.update_epsilon_sigma();
                //printf("epssum = %20.15f, sigmae = %20.15f\n", phen.get_epsilon_sum(), phen.get_sigmae());
            }
            phen.set_mu(phen.sample_norm_rng());
            //printf("new mu = %20.15f\n", phen.get_mu());
            phen.offset_epsilon(-phen.get_mu());
            //**BUG in original phen.epsilon_stats();

            // Shuffling of the markers on its own PRNG (see README/wiki)
            if (opt.shuffle_markers())
                phen.shuffle_midx(opt.mimic_hydra());

            phen.reset_m0();
            phen.reset_cass();
        }
        fflush(stdout);

        // Covariates
        if(C > 0){
            for(int covi = 0; covi < C; covi++){
                for (auto& phen : pmgr.get_phens()) {
                    // Update epsilon with respect to previous covariate effect
                    double delta = phen.get_cov_delta(covi);
                    phen.epsilon_update_cov(covi, delta);

                    double cov_num = phen.cov_dot_product(covi);
                    double cov_denom = phen.get_cov_denom(covi);
                    double delta_new = 0.0;

                    if(cov_denom > 0){
                        // Sample new covariate effect delta
                        delta_new = phen.sample_norm_rng(cov_num / cov_denom, 1.0 / cov_denom);
                    }

                    phen.set_cov_delta(covi, delta_new);
                    //printf("New delta[%d] = %0.6f\n", covi, delta_new);
                    //fflush(stdout);

                    // Update epsilon with respect to the new covariate effect
                    phen.epsilon_update_cov(covi, -delta_new);
                }
            }
        }
        //fflush(stdout);

        double dbetas[NPHEN * 3]; // [ dbeta:mave:msig | ... | ]

        double t_it_sync = 0.0;

        for (int mrki=0; mrki<Mm; mrki++) {

            bool share_mrk = false;
            for (int i=0; i<NPHEN*3; i++)
                dbetas[i] = 0.0;
            int mloc = 0;

            if (mrki < M) {

                mloc = pmgr.get_phens()[0].get_marker_local_index(mrki);
                const int mglo = S + mloc;
                const int mgrp = get_marker_group(mglo);
                //std::cout << "mloc = " << mloc << ", mglo = " << mglo << ", mgrp = " << mgrp << std::endl;

                int pheni = -1;

                for (auto& phen : pmgr.get_phens()) {

                    pheni += 1;

                    // Adav
                    if (phen.get_sigmag_for_group(mgrp) == 0.0) {
                        phen.set_marker_acum(mloc, 1.0);
                        phen.set_marker_beta(mloc, 0.0);
                        continue;
                    }

                    double beta   = phen.get_marker_beta(mloc);
                    double sige_g = phen.get_sigmae() / phen.get_sigmag_for_group(get_marker_group(mglo));
                    double sigg_e = 1.0 / sige_g;
                    double inv2sige = 1.0 / (2.0 * phen.get_sigmae());
                    //if (mrki < 3)
                        //printf("mrk = %3d has beta = %20.15f; sige_g = %20.15f = %20.15f / %20.15f\n", mloc, phen.get_marker_beta(mloc), sige_g, phen.get_sigmae(), phen.get_sigmag());

                    std::vector<double> denom = phen.get_denom();
                    std::vector<double> muk   = phen.get_muk();
                    std::vector<double> logl  = phen.get_logl();

                    for (int i=1; i<=K-1; ++i) {
                        denom.at(i-1) = (double)(N - 1) + sige_g * cvai[mgrp][i];
                        //printf("it %d, rank %d, m %d: denom[%d] = %20.15f, cvai = %20.15f\n", it, rank, mloc, i-1, denom.at(i-1), cvai[mgrp][i]);
                    }
                    double ts_debug_dp = MPI_Wtime();
                    double num = dot_product(mloc, phen.get_epsilon(), phen.get_marker_ave(mloc), phen.get_marker_sig(mloc));
                    double te_debug_dp = MPI_Wtime();
                    
                    //printf("RANK %d, time dot product: %10.3f\n", rank, (te_debug_dp - ts_debug_dp));

                    //printf("num = %20.15f\n", num);
                    num += beta * double(phen.get_nonas() - 1);

                    //printf("i:%d r:%d m:%d: num = %.17g, %20.15f, %20.15f\n", it, rank, mloc, num, phen.get_marker_ave(mloc), phen.get_marker_sig(mloc));

                    for (int i=1; i<=K-1; ++i)
                        muk.at(i) = num / denom.at(i - 1);

                    for (int i=0; i<K; i++) {
                        logl[i] = log(phen.get_pi_est(mgrp, i));
                        if (i>0)
                            logl[i] += -0.5 * log(sigg_e * double(phen.get_nonas() - 1) * cva[mgrp][i] + 1.0) + muk[i] * num * inv2sige;
                        //printf("logl[%d] = %20.15f\n", i, logl[i]);
                    }

                    double prob = phen.sample_unif_rng();

                    bool zero_acum = false;
                    double tmp1 = 0.0;
                    for (int i=0; i<K; i++) {
                        if (abs(logl[i] - logl[0]) > 700.0)
                            zero_acum = true;
                        tmp1 += exp(logl[i] - logl[0]);
                    }
                    zero_acum ? tmp1 = 0.0 : tmp1 = 1.0 / tmp1;
                    phen.set_marker_acum(mloc, tmp1);
                    //printf("i:%d r:%d m:%d p:%d  num = %20.15f, acum = %20.15f, prob = %20.15f\n", it, rank, mloc, pheni, num, phen.get_marker_acum(mloc), prob);
                    //fflush(stdout);
                    double dbeta = phen.get_marker_beta(mloc);

                    for (int i=0; i<K; i++) {
                        if (prob <= phen.get_marker_acum(mloc) || i == K - 1) {
                            if (i == 0) {
                                phen.set_marker_beta(mloc, 0.0);
                                //printf("@0@ i:%4d r:%4d m:%4d: beta reset to 0.0\n", it, rank, mloc);
                            } else {
                                phen.set_marker_beta(mloc, phen.sample_norm_rng(muk[i], phen.get_sigmae() / denom[i-1]));
                                //printf("@B@ i:%4d r:%4d m:%4d:  dbetat = %20.15f, muk[%4d] = %15.10f with prob=%15.10f <= acum = %15.10f, denom = %15.10f, sigmaE = %15.10f: beta = %15.10f\n", it, rank, mloc, phen.get_marker_beta(mloc) - dbeta, i, muk[i], prob, phen.get_marker_acum(mloc), denom[i-1], phen.get_sigmae(), phen.get_marker_beta(mloc));
                                //fflush(stdout);
                            }
                            phen.increment_cass(mgrp, i, 1);
                            //std::cout << "cass " << mgrp << " " << i << " = " << phen.get_cass_for_group(mgrp,i) << std::endl;
                            phen.set_comp(mloc, i);
                            break;
                        } else {
                            bool zero_inc = false;
                            for (int j=i+1; j<K; j++) {
                                if (abs(logl[j] - logl[i+1]) > 700.0)
                                    zero_inc = true;
                            }
                            if (!zero_inc) {
                                double esum = 0.0;
                                for (int k=0; k<logl.size(); k++)
                                    esum += exp(logl[k] - logl[i+1]);
                                phen.set_marker_acum(mloc, phen.get_marker_acum(mloc) + 1.0 / esum);
                            }
                        }
                    }

                    dbeta -= phen.get_marker_beta(mloc);

                    //printf("iteration %3d, rank %3d, marker %5d: dbeta = %20.15f, pheni = %d, %20.15f %20.15f\n", it, rank, mloc, dbeta, pheni, dbetas[pheni + 1], dbetas[pheni + 2]);

                    if (abs(dbeta) > 0.0) {
                        share_mrk = true;
                        dbetas[pheni * 3 + 0] = dbeta;
                        dbetas[pheni * 3 + 1] = phen.get_marker_ave(mloc);
                        dbetas[pheni * 3 + 2] = phen.get_marker_sig(mloc);
                        //printf("@¢@ i:%3d r:%3d m:%5d: dbeta = %20.15f, pheni = %d, %20.15f %20.15f\n", it, rank, mloc, dbeta, pheni, dbetas[pheni * 3 + 1], dbetas[pheni * 3 + 2]);
                    } else {
                        //printf("-¢- i:%3d r:%3d m:%5d: dbeta = %20.15f, pheni = %d, %20.15f %20.15f\n", it, rank, mloc, dbeta, pheni, dbetas[pheni * 3 + 1], dbetas[pheni * 3 + 2]);
                    }
                }
            }

            MPI_Barrier(MPI_COMM_WORLD);
            double ts_sync = MPI_Wtime();

            // Collect information on which tasks need to share its just processed marker;
            // Tasks with M < Mm, share_mrk is false by default.
            MPI_Allgather(&share_mrk,  1, MPI_C_BOOL,
                          recv_update, 1, MPI_C_BOOL,
                          MPI_COMM_WORLD);

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

            int totbytes = 0;

            int dis_bet[nranks], cnt_bet[nranks];
            int dis_meth[nranks], cnt_meth[nranks];
            int disp_bet = 0, disp_meth = 0;
            for (int i=0; i<nranks; i++) {
                if (recv_update[i]) {
                    cnt_bet[i] = NPHEN * 3;
                    cnt_meth[i] = N;
                    totbytes += mbytes;
                } else {
                    cnt_bet[i] = 0;
                    cnt_meth[i] = 0;
                }
                dis_bet[i] = disp_bet;
                dis_meth[i] = disp_meth;
                //printf("mrki = %d, recv_update[%d] = %s %d:%d\n", mrki, i, recv_update[i] ? "T" : "F", dis_bet[i], cnt_bet[i]);
                disp_bet += cnt_bet[i];
                disp_meth += cnt_meth[i];
            }

            //printf("rank %d bytes vs %d\n", totbytes, disp_meth);
            double ts_debug_rdb = MPI_Wtime();
            double recv_dbetas[disp_bet];
            MPI_Allgatherv(&dbetas, share_mrk ? NPHEN * 3 : 0, MPI_DOUBLE,
                           recv_dbetas, cnt_bet, dis_bet, MPI_DOUBLE, MPI_COMM_WORLD);

            double te_debug_rdb = MPI_Wtime();
            //printf("recv_dbetas gathered...\n");
            //fflush(stdout);

            double ts_debug_rm = MPI_Wtime();
            double* recv_meth = (double*) _mm_malloc(disp_meth * sizeof(double), 32);
            //printf("recv_meth allocated...\n");
            //fflush(stdout);
            check_malloc(recv_meth, __LINE__, __FILE__);
            //printf("recv_meth checked...\n");
            //fflush(stdout);
            size_t methix = size_t(mloc) * size_t(N);

            MPI_Allgatherv(&meth_data[methix], share_mrk ? size_t(N) : 0, MPI_DOUBLE,
                           recv_meth, cnt_meth, dis_meth, MPI_DOUBLE, MPI_COMM_WORLD);
            double te_debug_rm = MPI_Wtime();
            //printf("recv_meth gathered...\n");
            //fflush(stdout);
            double ts_debug_ue = MPI_Wtime();
            update_epsilon(cnt_bet, recv_dbetas, recv_meth);
            double te_debug_ue = MPI_Wtime();

            //if(rank == 0){
            //    printf("INFO: Time recv dbetas: %10.3f, Time recv meth: %10.3f, Time epsilon update: %10.3f\n", (te_debug_rdb - ts_debug_rdb), (te_debug_rm - ts_debug_rm), (te_debug_ue - ts_debug_ue));
            //}
            //printf("epsilon updated...\n");
            //fflush(stdout);
            MPI_Barrier(MPI_COMM_WORLD);
            double te_sync = MPI_Wtime();
            t_it_sync += te_sync - ts_sync;
            //printf("time %f \n", t_it_sync);
            //fflush(stdout);
            _mm_free(recv_meth);
            //printf("mrk %d end loop\n", mloc);
            //printf("-----------------------\n");
            fflush(stdout);

        } // End marker loop

        if(rank == 0){
            printf("INFO: sinc time: %10.3f\n", t_it_sync);
        }

        //exit(0);
        //continue;

        //MPI_Barrier(MPI_COMM_WORLD);
        double ts_debug_gp = MPI_Wtime();
        int pheni = -1;
        for (auto& phen : pmgr.get_phens()) {
            pheni++;

            phen.reset_beta_sqn_to_zero();
            for (int i=0; i<M; i++) {
                phen.increment_beta_sqn(get_marker_group(S + i), phen.get_marker_beta(i) * phen.get_marker_beta(i));
            }

            //for (int i=0; i<G; i++)
            //    printf("i:%d r:%d p:%d g:%d: beta_sqn[%d] = %20.15f\n", it, rank, pheni, i, i, phen.get_beta_sqn_for_group(i));

            double* beta_sqn_sum = (double*) malloc(G * sizeof(double));
            check_mpi(MPI_Allreduce(phen.get_beta_sqn()->data(),
                                    beta_sqn_sum,
                                    G,
                                    MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD), __LINE__, __FILE__);
            phen.set_beta_sqn(beta_sqn_sum);
            free(beta_sqn_sum);

            int* cass_sum = (int*) malloc(G * K * sizeof(int));
            check_mpi(MPI_Allreduce(phen.get_cass(),
                                    cass_sum,
                                    G * K,
                                    MPI_INT,
                                    MPI_SUM,
                                    MPI_COMM_WORLD), __LINE__, __FILE__);
            phen.set_cass(cass_sum);
            free(cass_sum);

            // Update global parameters
            //
            for (int i=0; i<G; i++) {

                // Skip empty groups
                if (mtotgrp.at(i) == 0)
                    continue;

                //if(rank == 0)
                //    printf(" !!! i:%d r:%d p:%d g:%d:  OLD sigmag = %20.15f\n", it, rank, pheni, i, phen.get_sigmag_for_group(i));

                //if (rank == 0)
                //    printf("i:%d r:%d p:%d g:%d:  m0 = %d - %d = %d\n", it, rank, pheni, i, mtotgrp.at(i), phen.get_cass_for_group(i, 0), mtotgrp.at(i) - phen.get_cass_for_group(i, 0));

                phen.set_m0_for_group(i, mtotgrp.at(i) - phen.get_cass_for_group(i, 0));

                // Skip groups with m0 being null or empty cass (adaV in action)
                if (phen.get_m0_for_group(i) == 0 || phen.get_cass_sum_for_group(i) == 0) {
                    phen.set_sigmag_for_group(i, 0.0);
                    continue;
                }
                //if(rank == 0)
                //    printf("m0g[%d] = %f, betasgn = %f\n", i, phen.get_m0_for_group(i), phen.get_beta_sqn_for_group(i));
                phen.set_sigmag_for_group(i, phen.sample_inv_scaled_chisq_rng(V0G + (double) phen.get_m0_for_group(i), (phen.get_beta_sqn_for_group(i) * (double) phen.get_m0_for_group(i) + V0G * S02G) / (V0G + (double) phen.get_m0_for_group(i))));
                //if(rank == 0)
                //    printf(" !!! i:%d r:%d p:%d g:%d:  NEW sigmag = %20.15f %20.15f\n", it, rank, pheni, i, phen.get_sigmag_for_group(i), phen.get_beta_sqn_for_group(i));

                phen.update_pi_est_dirichlet(i);
            }

            //continue;

            //if (rank == 0)
            //    phen.print_cass(mtotgrp);


            // Broadcast sigmaG of rank 0
            check_mpi(MPI_Bcast(phen.get_sigmag()->data(), phen.get_sigmag()->size(), MPI_DOUBLE, 0, MPI_COMM_WORLD), __LINE__, __FILE__);
            //for (int i=0; i<opt.get_ngroups(); i++) {
            //    printf("i:%d r:%d p:%d g:%d:  sigmag = %20.15f\n", it, rank, pheni, i, phen.get_sigmag_for_group(i));
            //}

            double e_sqn = phen.epsilon_sumsqr();
            //printf("i:%d r:%d p:%d  e_sqn = %20.15f\n", it, rank, pheni, e_sqn);

            //EO: sample sigmaE and broadcast the one from rank 0 to all the others
            phen.set_sigmae(phen.sample_inv_scaled_chisq_rng(V0E + (double)N, (e_sqn + V0E * S02E) / (V0E + (double)N)));

            double sigmae_r0 = phen.get_sigmae();
            check_mpi(MPI_Bcast(&sigmae_r0, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD), __LINE__, __FILE__);
            phen.set_sigmae(sigmae_r0);
            if (rank % 10 == 0) {
                printf("RESULT : i:%d r:%d p:%d  sum sigmaG = %20.15f  sigmaE = %20.15f\n", it, rank, pheni, phen.get_sigmag_sum(), phen.get_sigmae());
                //for (int i=0; i<opt.get_ngroups(); i++) {
                //    printf("i:%d r:%d p:%d  sigmaG = %20.15f  sigmaE = %20.15f\n", it, rank, pheni, phen.get_sigmag_for_group(i), phen.get_sigmae());
                //}
            }

            // Broadcast pi_est from rank 0
            for (int i=0; i<G; i++) {
                check_mpi(MPI_Bcast(phen.get_pi_est()->at(i).data(), phen.get_pi_est()->at(i).size(), MPI_DOUBLE, 0, MPI_COMM_WORLD), __LINE__, __FILE__);
            }
        }

        double te_it = MPI_Wtime();
        //if (rank == 0)
        //    printf("INFO : It %d  Global parameters time = %7.3f\n", it, te_it - ts_debug_gp);
        if (rank == 0)
            printf("RESULT : It %d  total proc time = %7.3f sec, with sync time = %7.3f\n", it, te_it - ts_it, t_it_sync);

        // Write output files
        if (it % opt.get_output_thin_rate() == 0) {
            const unsigned nthinned = it / opt.get_output_thin_rate() - 1;
            for (auto& phen : pmgr.get_phens()) {
                if (rank == 0) {
                    //phen.print_pi_est();
                    write_ofile_csv(*(phen.get_outcsv_fh()), it,  phen.get_sigmag(), phen.get_sigmae(), phen.get_m0_sum(), nthinned, phen.get_pi_est());
                    if(C > 0)
                        write_ofile_cov(*(phen.get_outcov_fh()), it,  phen.get_deltas(), nthinned);
                }
                write_ofile_h1(*(phen.get_outbet_fh()), rank, Mt, it, nthinned, S, M, phen.get_betas().data(), MPI_DOUBLE);
                write_ofile_h1(*(phen.get_outcpn_fh()), rank, Mt, it, nthinned, S, M, phen.get_comp().data(),  MPI_INTEGER);
            }
        }

    } // End iteration loop

    MPI_Barrier(MPI_COMM_WORLD);

    for (auto& phen : pmgr.get_phens())
        phen.close_output_files();

    double te_overall = MPI_Wtime();
    if (rank == 0)
        printf("OVERALL : Overall proc time = %7.3f sec\n",te_overall - ts_overall);
}

// counts[rank]: holds either NPHEN * 3 or 0

void Bayes::update_epsilon(const int* counts, const double* dbetas, const double* meth) {

    const int NPHEN =  pmgr.get_phens().size();

    int cnt_tot = 0;
    int methi = 0;
    for (int i=0; i<nranks; i++) {
        int cnt = counts[i];
        //printf("r:%d  cnt = %d\n", i, cnt);
        //fflush(stdout);
        if (cnt == 0) continue;
        assert(cnt == NPHEN * 3);
        //printf(" cnt > 0: task %d, get for rank %d cnt = %d dbetas:\n", rank, i, cnt);
        //fflush(stdout);
        int pheni = 0;
        for (auto& phen : pmgr.get_phens()) {
            //printf("r:%d p:%d - dbeta = %20.15f\n", rank, pheni, dbetas[cnt_tot + pheni * 3]);
            //fflush(stdout);
            size_t methix = size_t(methi) * size_t(N);
            if (dbetas[cnt_tot + pheni * 3] != 0.0)
                phen.update_epsilon(&dbetas[cnt_tot + pheni * 3], &meth[methix]);
            //phen.update_epsilon_sum();
            //printf(" - r:%d p:%d   epssum = %.17g, cnt_tot = %d\n", i, pheni, phen.get_epsilon_sum(), cnt_tot);
            //fflush(stdout);
            pheni += 1;
        }
        cnt_tot += cnt;
        methi    += 1;
    }
    fflush(stdout);
}

// EO, review!!
double Bayes::dot_product(const int mrki, double* __restrict__ phen, const double mu, const double sigma_inv) {
    size_t methix = size_t(mrki) * size_t(N);
    double* meth = &meth_data[methix];
    double dpa = 0.0;

#ifdef _OPENMP
#pragma omp parallel for simd schedule(static) reduction(+:dpa)
#endif
    for (int i=0; i<N; i++)
        dpa += (meth[i] - mu) * phen[i];

    return sigma_inv * dpa;
}


// Setup processing: load input files and define MPI task workload
void Bayes::setup_processing() {

    mbytes = N * sizeof(double);//(N %  4) ? (size_t) N /  4 + 1 : (size_t) N /  4;
    //load_genotype();
    load_methylation_data();

    int C = opt.get_cov_num();
    if(C > 0)
        load_covariates();

    pmgr.read_phen_files(opt, get_N(), get_M());

    check_processing_setup();

    if (rank == 0)
        printf("INFO   : output directory: %s\n", opt.get_out_dir().c_str());

    MPI_Barrier(MPI_COMM_WORLD);
    double ts = MPI_Wtime();
    pmgr.compute_markers_statistics(meth_data, get_N(), get_M(), mbytes);
    MPI_Barrier(MPI_COMM_WORLD);
    double te = MPI_Wtime();
    if (rank == 0)
        printf("INFO   : Time to compute the markers' statistics: %.2f seconds.\n", te - ts);

    if (opt.predict()) return;

    for (auto& phen : pmgr.get_phens()) {
        phen.set_Z(Z); // set covariates 
        phen.set_prng_m((unsigned int)(opt.get_seed() + (rank + 0)));
        if (opt.mimic_hydra()) {
            phen.set_prng_d((unsigned int)(opt.get_seed() + (rank + 0) * 1000));
        } else {
            phen.set_prng_d((unsigned int)(opt.get_seed() + (rank + 1) * 1000));
        }
    }

    read_group_index_file(opt.get_group_index_file());

    for (int i=0; i<Mt; i++) {
        mtotgrp.at(get_marker_group(i)) += 1;
    }
    //for (int i=0; i<G; i++)
    //    printf("mtotgrp at %d = %d\n", i, mtotgrp.at(i));
    //pmgr.display_markers_statistics(get_M());
}


void Bayes::check_openmp() {
#ifdef _OPENMP
#pragma omp parallel
    {
        if (rank == 0) {
            int nt = omp_get_num_threads();
            if (omp_get_thread_num() == 0)
                printf("INFO   : OMP parallel regions will use %d thread(s)\n", nt);
        }
    }
#else
    printf("WARNING: no OpenMP support!\n");
#endif
}

void Bayes::read_group_index_file(const std::string& file) {

    std::ifstream infile(file.c_str());
    if (! infile)
        throw ("Error: can not open the group file [" + file + "] to read. Use the --group-index-file option!");

    if (rank == 0)
        std::cout << "INFO   : Reading groups from " + file + "." << std::endl;

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

void Bayes::load_methylation_data(){

    double ts = MPI_Wtime();

    MPI_File methfh;
    const std::string methfp = opt.get_bed_file();

    check_mpi(MPI_File_open(MPI_COMM_WORLD, methfp.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &methfh),  __LINE__, __FILE__);

    const size_t size_bytes = size_t(M) * size_t(N) * sizeof(double);

    meth_data = (double*)_mm_malloc(size_bytes, 32);

    printf("INFO  : rank %d has allocated %zu bytes (%.3f GB) for raw data.\n", rank, size_bytes, double(size_bytes) / 1.0E9);

    // Offset to section of bed file to be processed by task
    MPI_Offset offset = size_t(0) + size_t(S) * size_t(N) * sizeof(double);

    // Gather the sizes to determine common number of reads
    size_t max_size_bytes = 0;
    check_mpi(MPI_Allreduce(&size_bytes, &max_size_bytes, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD), __LINE__, __FILE__);

    const int NREADS = size_t( ceil(double(max_size_bytes)/double(INT_MAX/2)) );
    size_t bytes = 0;
    mpi_file_read_at_all <double*> (size_t(M) * size_t(N), offset, methfh, MPI_DOUBLE, NREADS, meth_data, bytes);

    MPI_File_close(&methfh);
    
    double te = MPI_Wtime();

    MPI_Barrier(MPI_COMM_WORLD);

    if (rank >= 0)
        printf("INFO   : time to load methylation data = %.2f seconds.\n", te - ts);

}

// values should be separate with space delimiter
void Bayes::load_covariates(){ 

    std::string covfp = opt.get_cov_file();
    int C = opt.get_cov_num();
    if (C==0)
        return;

    double ts = MPI_Wtime();

    std::ifstream covf(covfp);
    std::string line; 
    std::regex re("\\s+");

    int line_i = 0;

    while (std::getline(covf, line)) // read the current line
    {
        if(line_i == 0) {
            line_i++;
            continue; //skip header
        } 

        int Cobs = 0;
        std::vector<double> entries;
        std::sregex_token_iterator iter(line.begin(), line.end(), re, -1);
        std::sregex_token_iterator re_end;
        
        ++iter; // skip individual ID
        ++iter; // skip family ID
        for ( ; iter != re_end; ++iter){
            entries.push_back(std::stod(*iter));
            Cobs++;
        }

        if (Cobs != C){
            printf("FATAL   : number of covariates = %d does not match to the specified number of covariates = %d.\n", Cobs, C);
            exit(EXIT_FAILURE);
        }   
        Z.push_back(entries);
        line_i++;       
    }

    if (rank == 0)
        printf("INFO   : Number of loaded covariate lines from .cov file = %d \n", line_i - 1);

    double te = MPI_Wtime();

    if (rank == 0)
        printf("INFO   : time to load covariates = %.2f seconds.\n", te - ts);
    
    // Normalize covariates
    for(int covi = 0; covi < C; covi++){
            
        long double cavg = 0.0;
        long double csig = 0.0;

        for (int i = 0; i < N; i++) {
            cavg += Z[i][covi];    
        }
        cavg = cavg / double(N);

        for (int i = 0; i < N; i++) {
            csig += ((Z[i][covi] - cavg) * (Z[i][covi] - cavg));            
        }
        csig = sqrt(csig / double(N));

        for (int i = 0; i < N; i++) {
            if(csig < 0.00000001)
                Z[i][covi] = 0;
            else 
                Z[i][covi] = (Z[i][covi] - cavg) / csig;
        }
    }
}

void Bayes::set_block_of_markers() {

    const int modu = Mt % nranks;
    const int size = Mt / nranks;

    Mm = Mt % nranks != 0 ? size + 1 : size;

    int len[nranks], start[nranks];
    int cum = 0;
    for (int i=0; i<nranks; i++) {
        len[i]  = i < modu ? size + 1 : size;
        start[i] = cum;
        cum += len[i];
    }
    //printf("cum %d vs %d  Mm = %d\n", cum, Mt, Mm);
    assert(cum == Mt);

    M = len[rank];
    S = start[rank];

    printf("INFO   : rank %4d has %d markers over tot Mt = %d, max Mm = %d, starting at S = %d\n", rank, M, Mt, Mm, S);
    //@todo: mpi check sum over tasks == Mt
}

void Bayes::print_cva() {
    printf("INFO   : mixtures for all groups:\n");
    for (int i=0; i<G; i++) {
        printf("         grp %2d: ", i);
        for (int j=0; j<K; j++) {
            printf("%7.5f ", cva[i][j]);
        }
        printf("\n");
    }
}

void Bayes::print_cvai() {
    printf("INFO   : inverse mixtures for all groups:\n");
    for (int i=0; i<G; i++) {
        printf("         grp %2d: ", i);
        for (int j=0; j<K; j++) {
            printf("%10.3f ", cvai[i][j]);
        }
        printf("\n");
    }
}
