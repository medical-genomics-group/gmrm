#include <iostream>
#include <mpi.h>
#include "options.hpp"
#include "bayes.hpp"
#include "phenotype.hpp"
#include "genotype.hpp"


int main(int argc, char *argv[]) {

    Options opt(argc, argv);

    BayesRR brr(opt);
    
    Genotype geno;


    MPI_Init(NULL, NULL);

    for (int iter = 1; iter <= 10; iter++) {
        
    }

    MPI_Finalize();

    return 0;
}
