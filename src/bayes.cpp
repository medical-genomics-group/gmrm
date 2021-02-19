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


void Bayes::process() {

    check_openmp();
    
    for (auto& phen : pmgr.get_phens()) {
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

    const int NPHEN =  pmgr.get_phens().size();
    const int UPDR  = 3; // Update every 3 markers for now; --option to come...
    MPI_Request req_update[UPDR];
    MPI_Request req_dbetas[UPDR];
    MPI_Status  req_status[UPDR];
    bool send_update[nranks];
    bool recv_update[nranks]; //* UPDR];


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
                    //    printf("mrk = %3d has beta = %20.15f; sige_g = %20.15f = %20.15f / %20.15f\n", mloc, phen.get_marker_beta(mloc), sige_g, phen.get_sigmae(), phen.get_sigmag());

                    std::vector<double> denom = phen.get_denom();
                    std::vector<double> muk   = phen.get_muk();
                    std::vector<double> logl  = phen.get_logl();
                    std::vector<int> comp     = phen.get_comp();
                    std::vector<std::vector<double>> pi_est = phen.get_pi_est();

                    for (int i=1; i<=K-1; ++i) {
                        denom.at(i-1) = (double)(N - 1) + sige_g * cvai[mgrp][i];
                        //printf("it %d, rank %d, m %d: denom[%d] = %20.15f, cvai = %20.15f\n", it, rank, mloc, i-1, denom.at(i-1), cvai[mgrp][i]);
                    }

                    double num = dot_product(mloc, phen.get_epsilon(), phen.get_marker_ave(mloc), phen.get_marker_sig(mloc));
                    
                    //printf("num = %20.15f\n", num);
                    num += beta * double(phen.get_nonas() - 1);
                    
                    //printf("i:%d r:%d m:%d: num = %.17g, %20.15f, %20.15f\n", it, rank, mloc, num, phen.get_marker_ave(mloc), phen.get_marker_sig(mloc));

                    for (int i=1; i<=K-1; ++i)
                        muk.at(i) = num / denom.at(i - 1);

                    for (int i=0; i<K; i++) {
                        logl[i] = log(pi_est[mgrp][i]);
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
                            comp[mloc] = i;
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
                    dbetas[pheni * 3 + 0] = dbeta;
                    dbetas[pheni * 3 + 1] = phen.get_marker_ave(mloc);
                    dbetas[pheni * 3 + 2] = phen.get_marker_sig(mloc);
                    //printf("iteration %3d, rank %3d, marker %5d: dbeta = %20.15f, pheni = %d, %20.15f %20.15f\n", it, rank, mloc, dbeta, pheni, dbetas[pheni + 1], dbetas[pheni + 2]);

                    if (abs(dbeta) > 0.0) {
                        share_mrk = true;
                        //printf("@¢@ i:%3d r:%3d m:%5d: dbeta = %20.15f, pheni = %d, %20.15f %20.15f\n", it, rank, mloc, dbeta, pheni, dbetas[pheni * 3 + 1], dbetas[pheni * 3 + 2]);
                    } else {
                        //printf("-¢- i:%3d r:%3d m:%5d: dbeta = %20.15f, pheni = %d, %20.15f %20.15f\n", it, rank, mloc, dbeta, pheni, dbetas[pheni * 3 + 1], dbetas[pheni * 3 + 2]);
                    }
                }
            }

            //MPI_Barrier(MPI_COMM_WORLD);
            double ts_sync = MPI_Wtime();

            // Collect information on which tasks need to share a marker;
            // Tasks with M < Mm, share_mrk is false by default.
            MPI_Allgather(&share_mrk,  1, MPI_C_BOOL,
                          recv_update, 1, MPI_C_BOOL,
                          MPI_COMM_WORLD);
            
            int dis_bet[nranks], cnt_bet[nranks];
            int dis_bed[nranks], cnt_bed[nranks];
            int disp_bet = 0, disp_bed = 0;
            for (int i=0; i<nranks; i++) {
                if (recv_update[i]) {
                    cnt_bet[i] = NPHEN * 3;
                    cnt_bed[i] = mbytes;
                } else {
                    cnt_bet[i] = 0;
                    cnt_bed[i] = 0;
                }
                dis_bet[i] = disp_bet; 
                dis_bed[i] = disp_bed;
                //printf("mrki = %d, recv_update[%d] = %s %d:%d\n", mrki, i, recv_update[i] ? "T" : "F", dis_bet[i], cnt_bet[i]);
                disp_bet += cnt_bet[i];
                disp_bed += cnt_bed[i];
            }
            
            double recv_dbetas[disp_bet];
            MPI_Allgatherv(&dbetas, share_mrk ? NPHEN * 3 : 0, MPI_DOUBLE,
                           recv_dbetas, cnt_bet, dis_bet, MPI_DOUBLE, MPI_COMM_WORLD);
            
            unsigned char* recv_bed = (unsigned char*) _mm_malloc(disp_bed, 32);
            check_malloc(recv_bed, __LINE__, __FILE__);
            MPI_Allgatherv(&bed_data[mloc * mbytes], share_mrk ? mbytes : 0, MPI_UNSIGNED_CHAR,
                           recv_bed, cnt_bed, dis_bed, MPI_UNSIGNED_CHAR, MPI_COMM_WORLD);
            
            update_epsilon(cnt_bet, recv_dbetas, recv_bed);
            
            MPI_Barrier(MPI_COMM_WORLD);
            double te_sync = MPI_Wtime();
            t_it_sync += te_sync - ts_sync;
            
            _mm_free(recv_bed);
            
        } // End marker loop

        //exit(0);
        //continue;

        //MPI_Barrier(MPI_COMM_WORLD);

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

                //printf(" !!! i:%d r:%d p:%d g:%d:  OLD sigmag = %20.15f\n", it, rank, pheni, i, phen.get_sigmag_for_group(i));

                //if (rank == 0)
                //    printf("i:%d r:%d p:%d g:%d:  m0 = %d - %d = %d\n", it, rank, pheni, i, mtotgrp.at(i), phen.get_cass_for_group(i, 0), mtotgrp.at(i) - phen.get_cass_for_group(i, 0));

                phen.set_m0_for_group(i, mtotgrp.at(i) - phen.get_cass_for_group(i, 0));

                // Skip groups with m0 being null or empty cass (adaV in action)
                if (phen.get_m0_for_group(i) == 0 || phen.get_cass_sum_for_group(i) == 0) {
                    phen.set_sigmag_for_group(i, 0.0);
                    continue;
                }

                phen.set_sigmag_for_group(i, phen.sample_inv_scaled_chisq_rng(V0G + (double) phen.get_m0_for_group(i), (phen.get_beta_sqn_for_group(i) * (double) phen.get_m0_for_group(i) + V0G * S02G) / (V0G + (double) phen.get_m0_for_group(i))));
                //printf(" !!! i:%d r:%d p:%d g:%d:  NEW sigmag = %20.15f %20.15f\n", it, rank, pheni, i, phen.get_sigmag_for_group(i), phen.get_beta_sqn_for_group(i));

                phen.update_pi_est_dirichlet(i);
            }

            //continue;

            //if (rank == 0)
            //    phen.print_cass(mtotgrp);


            // Broadcast sigmaG of rank 0
            check_mpi(MPI_Bcast(phen.get_sigmag()->data(), phen.get_sigmag()->size(), MPI_DOUBLE, 0, MPI_COMM_WORLD), __LINE__, __FILE__);
            for (int i=0; i<opt.get_ngroups(); i++) {
                //printf("i:%d r:%d p:%d g:%d:  sigmag = %20.15f\n", it, rank, pheni, i, phen.get_sigmag_for_group(i));
            }

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
                check_mpi(MPI_Bcast(phen.get_pi_est().at(i).data(), phen.get_pi_est().at(i).size(), MPI_DOUBLE, 0, MPI_COMM_WORLD), __LINE__, __FILE__);
            }
        }
        
        double te_it = MPI_Wtime();
        if (rank == 0)
            printf("RESULT : It %d  total proc time = %7.3f sec, with sync time = %7.3f\n", it, te_it - ts_it, t_it_sync);

    } // End iteration loop
}

// counts[rank]: holds either NPHEN * 3 or 0

void Bayes::update_epsilon(const int* counts, const double* dbetas, const unsigned char* bed) {

    const int NPHEN =  pmgr.get_phens().size();

    int cnt_tot = 0;
    int bedi = 0;
    for (int i=0; i<nranks; i++) {
        int cnt = counts[i];
        //printf("r:%d  cnt = %d\n", i, cnt);
        if (cnt == 0) continue;        
        assert(cnt == NPHEN * 3);
        //printf(" cnt > 0: task %d, get for rank %d cnt = %d dbetas:\n", rank, i, cnt);
        int pheni = 0;
        for (auto& phen : pmgr.get_phens()) {
            //printf(" - %20.15f\n", dbetas[cnt_tot + pheni * 3]);
            if (dbetas[cnt_tot + pheni * 3] != 0.0)
                phen.update_epsilon(&dbetas[cnt_tot + pheni * 3], &bed[bedi * mbytes]);
            //phen.update_epsilon_sum();
            //printf(" - r:%d p:%d   epssum = %.17g, cnt_tot = %d\n", i, pheni, phen.get_epsilon_sum(), cnt_tot);
            pheni += 1;
        }
        cnt_tot += cnt;
        bedi    += 1;
    }
}

// EO, review!!
double Bayes::dot_product(const int mloc, double* __restrict__ phen, const double mu, const double sigma_inv) {

    unsigned char* bed = &bed_data[mloc * mbytes];

#ifdef MANVEC
    __m256d luta, lutb; //, lutna;
    __m512d lutab, p42;
    __m256d p4   = _mm256_set1_pd(0.0);
    __m256d suma = _mm256_set1_pd(0.0);
    __m256d sumb = _mm256_set1_pd(0.0);
    __m512d sum42 = _mm512_set1_pd(0.0);

#ifdef _OPENMP
    //#pragma omp parallel for schedule(static) reduction(addpd4:suma,sumb)
#pragma omp parallel for schedule(static) reduction(addpd8:sum42)
#endif                              
    for (int j=0; j<mbytes; j++) {
        //luta  = _mm256_load_pd(&dotp_lut_a[bed[j] * 4]);
        //lutb  = _mm256_load_pd(&dotp_lut_b[bed[j] * 4]);
        //luta  = _mm256_load_pd(&dotp_lut_ab[bed[j] * 8]);
        //lutb  = _mm256_load_pd(&dotp_lut_ab[bed[j] * 8 + 4]);
        lutab = _mm512_load_pd(&dotp_lut_ab[bed[j] * 8]);
        //p4    = _mm256_load_pd(&phen[j * 4]);
        p42 = _mm512_broadcast_f64x4(_mm256_load_pd(&phen[j * 4]));
        ////lutna = _mm256_load_pd(&na_lut[mask4[j] * 4]); // phen = 0.0 on NAs!
        //luta  = _mm256_mul_pd(luta, p4);
        //lutb  = _mm256_mul_pd(lutb, p4);
        p42 = _mm512_mul_pd(p42, lutab);
        //suma  = _mm256_add_pd(suma, luta);
        //sumb  = _mm256_add_pd(sumb, lutb);
        sum42 = _mm512_add_pd(sum42, p42);
    }

    //return sigma_inv * 
    //    (suma[0] + suma[1] + suma[2] + suma[3] - mu * (sumb[0] + sumb[1] + sumb[2] + sumb[3]));
    return sigma_inv * 
        (sum42[0] + sum42[1] + sum42[2] + sum42[3] - mu * (sum42[4] + sum42[5] + sum42[6] + sum42[7]));

#else

    double dpa = 0.0;
    double dpb = 0.0;

#ifdef _OPENMP
    //#pragma omp parallel for schedule(static) reduction(addpd4:suma,sumb)
#pragma omp parallel for schedule(static) reduction(+:dpa,dpb)
#endif                              
    for (int i=0; i<mbytes; i++) {
#ifdef _OPENMP
#pragma omp simd aligned(dotp_lut_a,dotp_lut_b,phen:32)
#endif
        for (int j=0; j<4; j++) {
            dpa += dotp_lut_a[bed[i] * 4 + j] * phen[i * 4 + j];
            dpb += dotp_lut_b[bed[i] * 4 + j] * phen[i * 4 + j];
        }
    }

    return sigma_inv * (dpa - mu * dpb);

#endif

}


// Setup processing: load input files and define MPI task workload
void Bayes::setup_processing() {

    mbytes = (N %  4) ? (size_t) N /  4 + 1 : (size_t) N /  4;

    load_genotype();
    pmgr.read_phen_files(opt, get_N(), get_M());
    check_processing_setup();
    MPI_Barrier(MPI_COMM_WORLD);
    double ts = MPI_Wtime();
    pmgr.compute_markers_statistics(bed_data, get_N(), get_M(), mbytes);
    MPI_Barrier(MPI_COMM_WORLD);
    double te = MPI_Wtime();
    if (rank == 0)
        printf("INFO   : Time to compute the markers' statistics: %.2f seconds.\n", te - ts);

    for (auto& phen : pmgr.get_phens()) {
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


void Bayes::load_genotype() {

    MPI_File bedfh;
    const std::string bedfp = opt.get_bed_file();
    check_mpi(MPI_File_open(MPI_COMM_WORLD, bedfp.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &bedfh),  __LINE__, __FILE__);

    const size_t size_bytes = size_t(M) * size_t(mbytes) * sizeof(unsigned char);

    bed_data = (unsigned char*)_mm_malloc(size_bytes, 64);
    check_malloc(bed_data, __LINE__, __FILE__);
    printf("rank %d allocation %zu bytes (%.3f GB) for the raw data. Marker start at S = %d\n", rank, size_bytes, double(size_bytes) / 1.0E9, S);

    // Offset to section of bed file to be processed by task
    MPI_Offset offset = size_t(3) + size_t(S) * size_t(mbytes) * sizeof(unsigned char);

    // Gather the sizes to determine common number of reads
    size_t max_size_bytes = 0;
    check_mpi(MPI_Allreduce(&size_bytes, &max_size_bytes, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD), __LINE__, __FILE__); 

    const int NREADS = check_int_overflow(size_t(ceil(double(max_size_bytes)/double(INT_MAX/2))), __LINE__, __FILE__);
    size_t bytes = 0;
    mpi_file_read_at_all <unsigned char*> (size_bytes, offset, bedfh, MPI_UNSIGNED_CHAR, NREADS, bed_data, bytes);
    //@@@MPI_Barrier(MPI_COMM_WORLD);

    check_mpi(MPI_File_close(&bedfh), __LINE__, __FILE__);
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
    
    std::cout << "rank " << rank << " has " << M << " markers over Mt = " << Mt << ", Mm = " << Mm << ", starting at S = " << S << std::endl;

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
