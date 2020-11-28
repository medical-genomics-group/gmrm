#include <iostream>
#include <limits.h>
#include <cmath>
#include "bayes.hpp"
#include "utilities.hpp"

void Bayes::process() {
    
    for (auto& phen : pmgr.get_phens()) {
        phen.sample_sigmag_beta_rng(1.0, 1.0);
        printf("sample sigmag = %20.15f\n", phen.get_sigmag());
    }

    for (unsigned int it = 1; it <= opt.get_iterations(); it++) {
        
        for (auto& phen : pmgr.get_phens()) {
            phen.set_midx();
            phen.offset_epsilon(phen.get_mu());
            phen.epsilon_stats();
            printf("epssum = %20.15f, sigmae = %20.15f\n", phen.get_epssum(), phen.get_sigmae());
            phen.sample_mu_norm_rng();
            printf("new mu = %20.15f\n", phen.get_mu());
            phen.offset_epsilon(-phen.get_mu());

            // Important that all phens to the shuffling of the
            // index array of markers midx to keep consistency of their prng

            if (opt.shuffle_markers())
                phen.shuffle_midx();
            for (int i=0; i<10; i++)
                printf("it %4d  midx[%7d] = %7d\n", it, i, phen.get_midx()[i]);

        }
        


    } // End iteration loop
}


// Setup processing: load input files and define MPI task workload
void Bayes::setup_processing() {

    mrk_bytes = (N %  4) ? (size_t) N /  4 + 1 : (size_t) N /  4;
    mrk_uints = (N % 16) ? (size_t) N / 16 + 1 : (size_t) N / 16;

    load_genotype();
    pmgr.read_phen_files(opt, get_N(), get_M());
    check_processing_setup();

    // Compute phenotype-dependent markers' statistics
    pmgr.compute_markers_statistics(bed_data, get_N(), get_M(), mrk_bytes);

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

    const size_t size_bytes = size_t(M) * size_t(mrk_bytes) * sizeof(unsigned char);

    bed_data = (unsigned char*)_mm_malloc(size_bytes, 64);
    check_malloc(bed_data, __LINE__, __FILE__);
    printf("rank %d allocation %zu bytes (%.3f GB) for the raw data.\n", rank, size_bytes, double(size_bytes) / 1.0E9);

    // Offset to section of bed file to be processed by task
    MPI_Offset offset = size_t(3) + size_t(S) * size_t(mrk_bytes) * sizeof(unsigned char);

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
    M = rank < modu ? size + 1 : size;
    std::cout << "rank " << rank << " has " << M << " markers over Mt = " << Mt << std::endl;

    //@todo: mpi check sum over tasks == Mt
}
