#include <iostream>
#include <fstream>
#include <regex>
#include <cmath>
#include <immintrin.h>
#include <mpi.h>
#include "utilities.hpp"
#include "phenotype.hpp"
#include "dotp_lut.hpp"
#include "na_lut.hpp"

//double _mm512_mask_reduce_add_pd (__mmask8 k, __m512d a)
//gen _mm256_load_pd 

// Compute mean and associated standard deviation for markers
// for each of the phenotypes (stats are NA dependent)
// ! one byte of bed  contains information for 4 individuals
// ! one byte of phen contains information for 8 individuals
void PhenMgr::compute_markers_statistics(const unsigned char* bed, const int N, const int M, const int mbytes) {

    for (auto& phen : get_phens()) {

        phen.print_info();
        const std::vector<unsigned char> mask4 = phen.get_mask4();
        double* mave = (double*) _mm_malloc(size_t(M) * sizeof(double), 64);  check_malloc(mave, __LINE__, __FILE__);
        double* msig = (double*) _mm_malloc(size_t(M) * sizeof(double), 64);  check_malloc(msig, __LINE__, __FILE__);

        std::cout << " - loop over " << M << " markers " << std::endl;
        double start = MPI_Wtime();
        for (int i=0; i<M; i++) {
            __m256d suma = _mm256_set1_pd(0.0);
            __m256d sumb = _mm256_set1_pd(0.0);
            __m256d luta, lutb, lutna;
            for (int j=0; j<mbytes; j++) {
                luta  = _mm256_load_pd(&dotp_lut_a[bed[i*mbytes + j] * 4]);
                lutb  = _mm256_load_pd(&dotp_lut_b[bed[i*mbytes + j] * 4]);
                lutna = _mm256_load_pd(&na_lut[mask4[j] * 4]);
                luta  = _mm256_mul_pd(luta, lutna);
                lutb  = _mm256_mul_pd(lutb, lutna);
                suma  = _mm256_add_pd(suma, luta);
                sumb  = _mm256_add_pd(sumb, lutb);
            }
            double asum = suma[0] + suma[1] + suma[2] + suma[3];
            double bsum = sumb[0] + sumb[1] + sumb[2] + sumb[3];
            double avg  = asum / bsum;

            __m256d vave = _mm256_set1_pd(-avg);
            __m256d sums = _mm256_set1_pd(0.0);
            for (int j=0; j<mbytes; j++) {
                luta  = _mm256_load_pd(&dotp_lut_a[bed[i*mbytes + j] * 4]);
                lutb  = _mm256_load_pd(&dotp_lut_b[bed[i*mbytes + j] * 4]);
                lutna = _mm256_load_pd(&na_lut[mask4[j] * 4]);
                luta  = _mm256_add_pd(luta, vave);    // - mu
                luta  = _mm256_mul_pd(luta, lutb);    // M -> 0.0
                luta  = _mm256_mul_pd(luta, lutna);   // NAs
                luta  = _mm256_mul_pd(luta, luta);    // ^2
                sums  = _mm256_add_pd(sums, luta);    // sum
            }
            double sig = sqrt((sums[0] + sums[1] + sums[2] + sums[3]) / (double(phen.get_nonas()) - 1.0));
            
            mave[i] = avg;
            msig[i] = sig;

            if (i < 10)
                printf("avg for marker %d = %20.15f +/- %20.15f %20.15f (%6.1f / %d) %20.15f\n", i, avg, sig, 1.0 / sig, asum, phen.get_nonas() - 1, sums[0] + sums[1] + sums[2] + sums[3]); 

        }
        double end = MPI_Wtime();
        std::cout << "The process took " << end - start << " seconds to run." << std::endl;
    }
}

void PhenMgr::print_info() {
    for (auto& phen : get_phens())
        phen.print_info();
}

void PhenMgr::read_phen_files(const Options& opt) {
    std::vector<std::string> phen_files = opt.get_phen_files();
    for (auto fp = phen_files.begin(); fp != phen_files.end(); ++fp) {
        if (opt.verbosity_level(3))
            std::cout << "Reading phenotype file: " << *fp << std::endl;
        Phenotype phenotype(*fp, opt);
        phens.push_back(phenotype);
    }
}


// Read phenotype file
// Assume PLINK format: Family ID, Individual ID, Phenotype
// One row per individual
void Phenotype::read_file(const Options& opt) {

    std::ifstream infile(filepath);
    std::string line;
    std::regex re("\\s+");

    if (infile.is_open()) {
        int line_n = 0;
        nonas = 0, nas = 0;
        while (getline(infile, line)) {            
            int m8 = line_n % 8;
            int m4 = line_n % 4;
            if (m8 == 0)  mask8.push_back(0xFF);
            if (m4 == 0)  mask4.push_back(0b00001111);

            std::sregex_token_iterator first{line.begin(), line.end(), re, -1}, last;
            std::vector<std::string> tokens{first, last};
            if (tokens[2] == "NA") {
                nas += 1;
                data.push_back(NAN);
                if (opt.verbosity_level(3)) {
                    std::cout << " ... found NA on line " << line_n << ", m8 = " << m8 << " on byte " << int(line_n / 8) << std::endl;
                    std::cout << " ... found NA on line " << line_n << ", m4 = " << m4 << " on byte " << int(line_n / 4) << std::endl;
                }
                mask8.at(int(line_n / 8)) &= ~(0b1 << m8);
                mask4.at(int(line_n / 4)) &= ~(0b1 << m4);
            } else {
                nonas += 1;
                data.push_back(atof(tokens[2].c_str()));
            }

            line_n += 1;

            if (opt.verbosity_level(3)) {
                // Print fully handled bytes
                if (line_n % 8 == 0 && line_n > 7 && line_n < 30) {
                    std::cout << "mask8[" << int(line_n / 8) - 1 << "] = " << unsigned(mask8.at(int(line_n / 8) - 1)) << std::endl;
                }
                if (line_n % 4 == 0 && line_n > 3 && line_n < 30) {
                    std::cout << "mask4[" << int(line_n / 4) - 1 << "] = " << unsigned(mask4.at(int(line_n / 4) - 1)) << std::endl;
                }
            }
        }
        infile.close();

        // Set last bits to 0 if ninds % 4 != 0
        const int m4 = line_n % 4;
        if (m4 != 0) {
            std::cout << "Setting last " << m4 << " bits to NAs" << std::endl;
            std::cout << "fatal: missing implementation" << std::endl;
            exit(1);
        }

    } else {
        std::cout << "FATAL: could not open phenotype file: " << filepath << std::endl;
        exit(EXIT_FAILURE);
    }


    //assert(nonas + nas == numInds);
    
}

void Phenotype::print_info() const {
    printf("INFO   : %s has %d NAs and %d non-NAs.\n", get_filepath().c_str(), nas, nonas);
}
