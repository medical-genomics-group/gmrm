#include <iostream>
//#include <mpi.h>
#include "options.hpp"
#include "phenotype.hpp"
#include "bayes.hpp"
#include "dimensions.hpp"

int main(int argc, char* argv[]) {

    //MPI_Init(NULL, NULL);
    const Options   opt(argc, argv);
    std::cout << "Options created."<<"\n";

    const Dimensions dims(opt);
    std::cout << "Dimensions created." << "\n";
    printf("There are %d samples and %d markers. \n", dims.get_nt(), dims.get_mt());

    BayesRR brr(opt, dims);
    std::cout << "BayesRR constructor created." << "\n";

    //if (opt.predict()) {
    //    brr.predict();
    //}
    //else {
    //    brr.process();
    //};
    brr.process();

    //MPI_Finalize();

    return 0;
}
