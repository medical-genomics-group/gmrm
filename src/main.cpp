#include <iostream>
#include <mpi.h>
#include "options.hpp"
#include "phenotype.hpp"
#include "bayes.hpp"
#include "dimensions.hpp"

int main(int argc, char *argv[]) {

    MPI_Init(NULL, NULL);

    const Options    opt(argc, argv);
    const Dimensions dims(opt);

    BayesRR brr(opt, dims);
    brr.process();
    
    std::cout << "__END_PROCESSING__" << std::endl;
    
    MPI_Finalize();
    return 0;
}
