#include <iostream>
#include <limits.h>
#include <cmath>
#include <immintrin.h>
#include "bayes.hpp"
#include "utilities.hpp"
#include "dotp_lut.hpp"


void Bayes::process() {
    
    for (auto& phen : pmgr.get_phens()) {
        phen.set_midx();
        phen.set_sigmag(phen.sample_beta_rng(1.0, 1.0));
        printf("sample sigmag = %20.15f\n", phen.get_sigmag());

        phen.set_pi_est(pi_prior);
    }

    const int NPHEN =  pmgr.get_phens().size();
    const int UPDR  = 3; // Update every 3 markers for now; --option to come...
    MPI_Request req_update[UPDR];
    MPI_Request req_dbetas[UPDR];
    MPI_Status  req_status[UPDR];
    bool send_update[nranks];
    bool recv_update[nranks]; //* UPDR];



    for (unsigned int it = 1; it <= opt.get_iterations(); it++) {
        
        int pidx = 0;
        for (auto& phen : pmgr.get_phens()) {
            pidx += 1;
            //printf("phen %d has mu = %20.15f\n", pidx, phen.get_mu());
            phen.offset_epsilon(phen.get_mu());
            phen.epsilon_stats();
            //printf("epssum = %20.15f, sigmae = %20.15f\n", phen.get_epssum(), phen.get_sigmae());
            phen.set_mu(phen.sample_norm_rng());
            //printf("new mu = %20.15f\n", phen.get_mu());
            phen.offset_epsilon(-phen.get_mu());
            //**BUG in original phen.epsilon_stats();
            
            // Important that all phens to the shuffling of the
            // index array of markers midx to keep consistency of their prng
            if (opt.shuffle_markers())
                phen.shuffle_midx();
            //for (int i=0; i<10; i++)
            //    printf("it %4d  midx[%7d] = %7d\n", it, i, phen.get_midx()[i]);
        }

        for (int mrki=0; mrki<Mm; mrki++) {

            if (mrki < M) {


                const int mloc = pmgr.get_phens()[0].get_marker_local_index(mrki);
                const int mglo = S + mloc;
                const int mgrp = get_marker_group(mglo);
                //std::cout << "mloc = " << mloc << ", mglo = " << mglo << ", mgrp = " << mgrp << std::endl;


                bool    share_mrk = false;
                double  dbetas[NPHEN];


                int pheni = -1;



                for (auto& phen : pmgr.get_phens()) {

                    pheni += 1;

                    double beta   = phen.get_marker_beta(mloc);
                    double sige_g = phen.get_sigmae() / phen.get_sigmag();
                    double sigg_e = 1.0 / sige_g;
                    double inv2sige = 1.0 / (2.0 * phen.get_sigmae());
                    //if (mrki < 3)
                    //    printf("mrk = %3d has beta = %20.15f; sige_g = %20.15f = %20.15f / %20.15f\n", mloc, phen.get_marker_beta(mloc), sige_g, phen.get_sigmae(), phen.get_sigmag());
                    
                    std::vector<double> denom = phen.get_denom();
                    std::vector<double> muk   = phen.get_muk();
                    std::vector<double> logl  = phen.get_logl();
                    std::vector<int> comp     = phen.get_comp();
                    std::vector<std::vector<double>> pi_est = phen.get_pi_est();
                    std::vector<std::vector<int>> cass = phen.get_cass();
                    

                    for (int i=1; i<=K-1; ++i) {
                        denom.at(i-1) = (double)(N - 1) + sige_g * cvai[mgrp][i];
                        //printf("it %d, rank %d, m %d: denom[%d] = %20.15f, cvai = %20.15f\n", it, rank, mloc, i-1, denom.at(i-1), cvai[mgrp][i]);
                    }
                    
                    double num = dot_product(mloc, phen.get_epsilon(), phen.get_marker_ave(mloc), phen.get_marker_sig(mloc));
                    //printf("num = %20.15f\n", num);
                    num += beta * double(phen.get_nonas() - 1);
                    //printf("it %d, rank %d, mark %d: num = %22.15f, %20.15f, %20.15f\n", it, rank, mloc, num, phen.get_marker_ave(mloc), 1.0 / phen.get_marker_sig(mloc));

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
                    //printf("it %d, rank %d, marker %d, acum = %20.15f, prob = %20.15f\n", it, rank, mloc, phen.get_marker_acum(mloc), prob);

                    double dbeta = phen.get_marker_beta(mloc);

                    for (int i=0; i<K; i++) {
                        if (prob <= phen.get_marker_acum(mloc) || i == K - 1) {
                            if (i == 0) {
                                phen.set_marker_beta(mloc, 0.0);
                            } else {
                                phen.set_marker_beta(mloc, phen.sample_norm_rng(muk[i], phen.get_sigmae() / denom[i-1]));
                                //printf("@B@ it=%4d, rank=%4d, mloc=%4d: muk[%4d] = %15.10f with prob=%15.10f <= acum = %15.10f, denom = %15.10f, sigmaE = %15.10f: beta = %15.10f\n", it, rank, mloc, i, muk[i], prob, phen.get_marker_acum(mloc), denom[i-1], phen.get_sigmae(), phen.get_marker_beta(mloc));
                                fflush(stdout);
                            }
                            cass[mgrp][i] += 1;
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
                    printf("iteration %3d, marker %5d: dbeta = %20.15f, pheni = %d\n", it, mloc, dbeta, pheni);
                    dbetas[pheni] = dbeta;
                    if (abs(dbeta) > 0.0) share_mrk = true;                    
                }


                // Collect information on which tasks need to share a marker
                MPI_Allgather(&share_mrk,  1, MPI_C_BOOL,
                              recv_update, 1, MPI_C_BOOL,
                              MPI_COMM_WORLD);

                int nupdate = 0;
                int disps[nranks], counts[nranks];
                int dis_bed[nranks], cnt_bed[nranks];
                int disp = 0, disp_bed = 0;
                for (int i=0; i<nranks; i++) {
                    if (recv_update[i]) {
                        counts[i]  = NPHEN;
                        cnt_bed[i] = mbytes;
                    } else {
                        counts[i]  = 0;
                        cnt_bed[i] = 0;
                    }
                    disps[i] = disp; 
                    dis_bed[i] = disp_bed;
                    printf("mrki = %d, recv_update[%d] = %s %d:%d\n", mrki, i, recv_update[i] ? "true" : "false", disps[i], counts[i]);
                    disp     += counts[i];
                    disp_bed += cnt_bed[i];
                }

                // collect dbetas
                double recv_dbetas[disp];
                MPI_Allgatherv(&dbetas, share_mrk ? NPHEN : 0, MPI_DOUBLE,
                               recv_dbetas, counts, disps, MPI_DOUBLE, MPI_COMM_WORLD);
                int cnt = 0;
                for (int i=0; i<nranks; i++) {
                    printf("from task %d:\n", i);
                    if (counts[i]) {
                        for (int j=0; j<counts[i]; j++) {
                            for (int k=0; k<NPHEN; k++) {
                                printf(" %d -> %20.15f\n", k, recv_dbetas[cnt]);
                                cnt++;
                            }
                        }
                    }
                    printf("\n");
                }
                fflush(stdout);

                // collect bed data
                unsigned char* recv_bed = (unsigned char*) _mm_malloc(disp * mbytes, 64);
                check_malloc(recv_bed, __LINE__, __FILE__);
                MPI_Allgatherv(&bed_data[mloc * mbytes], mbytes, MPI_UNSIGNED_CHAR,
                               recv_bed, cnt_bed, dis_bed, MPI_UNSIGNED_CHAR, MPI_COMM_WORLD);
                _mm_free(recv_bed);

            } else {
                std::cout << "FATAL  : handle me please!" << std::endl;
                exit(1);
            }
        }

    } // End iteration loop
}

double Bayes::dot_product(const int mloc, double* phen, const double mu, const double sigma) {

    unsigned char* bed = &bed_data[mloc * mbytes];

    __m256d luta, lutb, lutna;
    __m256d p4   = _mm256_set1_pd(0.0);
    __m256d suma = _mm256_set1_pd(0.0);
    __m256d sumb = _mm256_set1_pd(0.0);

    for (int j=0; j<mbytes; j++) {
        luta  = _mm256_load_pd(&dotp_lut_a[bed[j] * 4]);
        lutb  = _mm256_load_pd(&dotp_lut_b[bed[j] * 4]);
        p4    = _mm256_load_pd(&phen[j * 4]);
        //lutna = _mm256_load_pd(&na_lut[mask4[j] * 4]); // phen = 0.0 on NAs!
        luta  = _mm256_mul_pd(luta, p4);
        lutb  = _mm256_mul_pd(lutb, p4);
        suma  = _mm256_add_pd(suma, luta);
        sumb  = _mm256_add_pd(sumb, lutb);
    }

    return (1.0 / sigma) * 
        (suma[0] + suma[1] + suma[2] + suma[3] - mu * (sumb[0] + sumb[1] + sumb[2] + sumb[3]));
}


// Setup processing: load input files and define MPI task workload
void Bayes::setup_processing() {

    mbytes = (N %  4) ? (size_t) N /  4 + 1 : (size_t) N /  4;

    load_genotype();
    pmgr.read_phen_files(opt, get_N(), get_M());
    check_processing_setup();
    pmgr.compute_markers_statistics(bed_data, get_N(), get_M(), mbytes);

    // All phenotypes get the same seed, and will sollicitate their own PRN
    // generator exactly the same way, so that the PRNG state is consistent
    for (auto& phen : pmgr.get_phens()) {
        phen.set_rng((unsigned int)(opt.get_seed() + rank*1000));
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
    printf("rank %d allocation %zu bytes (%.3f GB) for the raw data.\n", rank, size_bytes, double(size_bytes) / 1.0E9);

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
    M  = rank < modu ? size + 1 : size;
    Mm = Mt % nranks != 0 ? size + 1 : size;
    std::cout << "rank " << rank << " has " << M << " markers over Mt = " << Mt << ", Mm = " << Mm << std::endl;

    //@todo: mpi check sum over tasks == Mt
}
