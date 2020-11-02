#pragma once
#include <assert.h>

void check_malloc(const void* ptr, const int linenumber, const char* filename);

void check_mpi(const int error, const int linenumber, const char* filename);

int check_int_overflow(const size_t n, const int linenumber, const char* filename);


// MPI_File_read_at_all handling count argument larger than INT_MAX
//
template <typename T>
void mpi_file_read_at_all(const size_t N, const MPI_Offset offset, const MPI_File fh, const MPI_Datatype MPI_DT, const int NREADS, T buffer, size_t &bytes) {

    if (NREADS == 0) return;

    int dtsize = 0;
    MPI_Type_size(MPI_DT, &dtsize);
    assert(dtsize == sizeof(buffer[0]));

    int SPLIT_ON = check_int_overflow(size_t(ceil(double(N)/double(NREADS))), __LINE__, __FILE__);
    int count = SPLIT_ON;
    bytes = 0;
        
    MPI_Status status;
    for (uint i=0; i<NREADS; ++i) {
            
        const size_t iim = size_t(i) * size_t(SPLIT_ON);

        // Last iteration takes only the leftover
        if (i == NREADS-1) count = check_int_overflow(N - iim, __LINE__, __FILE__);

        check_mpi(MPI_File_read_at(fh, offset + iim * size_t(dtsize), &buffer[iim], count, MPI_DT, &status), __LINE__, __FILE__);

        bytes += size_t(count) * size_t(dtsize);
    }
}


// MPI_File_write_at_all handling count argument larger than INT_MAX
//
template <typename T>
void mpi_file_write_at_all(const size_t N, MPI_Offset offset, MPI_File fh, MPI_Datatype MPI_DT, const int NWRITES, T buffer) {

    if (NWRITES == 0) return;

    int dtsize = 0;
    MPI_Type_size(MPI_DT, &dtsize);
    assert(dtsize == sizeof(buffer[0]));


    int SPLIT_ON = check_int_overflow(size_t(ceil(double(N)/double(NWRITES))), __LINE__, __FILE__);
    int count = SPLIT_ON;

    MPI_Status status;
    for (uint i=0; i<NWRITES; ++i) {
        
        const size_t iim = size_t(i) * size_t(SPLIT_ON);

        // Last iteration takes only the leftover
        if (i == NWRITES-1) count = check_int_overflow(N - iim, __LINE__, __FILE__);

        check_mpi(MPI_File_write_at_all(fh, offset + iim * size_t(dtsize), &buffer[iim], count, MPI_DT, &status), __LINE__, __FILE__);
    }
}
