#include <iostream>
#include <limits.h>
#include <cmath>
#include "bayes.hpp"
#include "utilities.hpp"


// Setup processing: load input files and define MPI task workload
void Bayes::setup_processing() {

    mrk_bytes = (Nt %  4) ? (size_t) Nt /  4 + 1 : (size_t) Nt /  4;
    mrk_uints = (Nt % 16) ? (size_t) Nt / 16 + 1 : (size_t) Nt / 16;

    set_block_of_markers();

    load_genotype();

}

void Bayes::load_genotype() {

    MPI_File bedfh;
    const std::string bedfp = opt.get_bed_file();
    check_mpi(MPI_File_open(MPI_COMM_WORLD, bedfp.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &bedfh),  __LINE__, __FILE__);

    const size_t size_bytes = size_t(M) * size_t(mrk_bytes) * sizeof(char);

    bed_data = (char*)_mm_malloc(size_bytes, 64);
    check_malloc(bed_data, __LINE__, __FILE__);
    printf("rank %d allocation %zu bytes (%.3f GB) for the raw data.\n", rank, size_bytes, double(size_bytes) / 1.0E9);

    // Offset to section of bed file to be processed by task
    MPI_Offset offset = size_t(3) + size_t(S) * size_t(mrk_bytes) * sizeof(char);

    // Gather the sizes to determine common number of reads
    size_t max_size_bytes = 0;
    check_mpi(MPI_Allreduce(&size_bytes, &max_size_bytes, 1, MPI_UNSIGNED_LONG_LONG, MPI_MAX, MPI_COMM_WORLD), __LINE__, __FILE__); 

    const int NREADS = check_int_overflow(size_t(ceil(double(max_size_bytes)/double(INT_MAX/2))), __LINE__, __FILE__);
    size_t bytes = 0;
    mpi_file_read_at_all <char*> (size_bytes, offset, bedfh, MPI_CHAR, NREADS, bed_data, bytes);

    MPI_Barrier(MPI_COMM_WORLD);

    check_mpi(MPI_File_close(&bedfh), __LINE__, __FILE__);
}

void Bayes::set_block_of_markers() {
    const int modu = Mt % nranks;
    const int size = Mt / nranks;
    M = rank < modu ? size + 1 : size;
    std::cout << "rank " << rank << " has " << M << " markers over Mt = " << Mt << std::endl;
    //@todo: mpi check sum over tasks == Mt
}

// Check for minimal setup: a bed file + a dim file + phen file(s)
void Bayes::check_options() {
    
    std::cout << "will check passed options for completeness" << std::endl;
    
    if (opt.get_bed_file() == "") {
        std::cout << "FATAL  : no bed file provided! Please use the --bedfile option." << std::endl;
        exit(EXIT_FAILURE);
    }
    std::cout << "  bed file: OK - " << opt.get_bed_file() << "\n";

    if (opt.get_dim_file() == "") {
        std::cout << "FATAL  : no dim file provided! Please use the --dimfile option." << std::endl;
        exit(EXIT_FAILURE);
    }
    std::cout << "  dim file: OK - " << opt.get_dim_file() << "\n";

    if (opt.count_phen_files() == 0) {
        std::cout << "FATAL  : no phen file(s) provided! Please use the --phenfile option." << std::endl;
        exit(EXIT_FAILURE);
    }
    std::cout << "  phen file(s): OK - " << opt.count_phen_files() << " files passed.\n";
    opt.list_phen_files();
}
