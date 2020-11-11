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
    const PhenMgr    pmgr(opt);

    BayesRR brr(opt, pmgr, dims);

    MPI_Finalize();
    return 0;
}
