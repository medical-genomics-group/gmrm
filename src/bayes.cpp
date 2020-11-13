#include <iostream>
#include <limits.h>
#include <cmath>
#include "bayes.hpp"
#include "utilities.hpp"
#include <boost/range/algorithm.hpp>

void Bayes::process() {

    for (int i=0; i<10; i++)
        printf("2: %6d\n", midx[i]);
    if (opt.shuffle_markers())  shuffle_markers();
    for (int i=0; i<10; i++)
        printf("1: %6d\n", midx[i]);
    if (opt.shuffle_markers())  shuffle_markers();
    for (int i=0; i<10; i++)
        printf("2: %6d\n", midx[i]);
}

void Bayes::shuffle_markers() {
    boost::uniform_int<> unii(0, M-1);
    boost::variate_generator< boost::mt19937&, boost::uniform_int<> > generator(dist.get_rng(), unii);
    boost::range::random_shuffle(midx, generator);
}

// Setup processing: load input files and define MPI task workload
void Bayes::setup_processing() {

    mrk_bytes = (Nt %  4) ? (size_t) Nt /  4 + 1 : (size_t) Nt /  4;
    mrk_uints = (Nt % 16) ? (size_t) Nt / 16 + 1 : (size_t) Nt / 16;

    set_block_of_markers();
    load_genotype();
    check_processing_setup();

    // Compute phenotype-dependent markers' statistics
    pmgr.compute_markers_statistics(bed_data, get_N(), get_M(), mrk_bytes);

    dist.set_rng((unsigned int)(opt.get_seed() + rank*1000));
}

void Bayes::check_processing_setup() {
    for (const auto& phen : pmgr.get_phens()) {
        const int Np = phen.get_nas() + phen.get_nonas();
        if (Np != Nt) {
            std::cout << "Fatal: Nt = " << Nt << " while phen file " << phen.get_filepath() << " has " << Np << " individuals!" << std::endl;
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

    // List of markers
    for (int i=0; i<M; ++i) midx.push_back(i);
}
