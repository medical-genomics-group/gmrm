#include <iostream>
#include <mpi.h>
#include "options.hpp"
#include "phenotype.hpp"
#include "bayes.hpp"
#include "dimensions.hpp"

int main(int argc, char *argv[]) {

    const Options opt(argc, argv);
    const PhenMgr phen_mgr(opt);

    MPI_Init(NULL, NULL);

    const Dimensions dims(opt);

    BayesRR brr(opt, phen_mgr, dims);

    for (int iter = 1; iter <= 10; iter++) {
        
    }

    MPI_Finalize();

    return 0;
}
