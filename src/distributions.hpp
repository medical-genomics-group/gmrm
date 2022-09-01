#pragma once
#include <random>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <vector>
#include <chrono>

#include <boost/random.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#define BOOST_UBLAS_NDEBUG 1

#include "pdflib.hpp"


using namespace std;



class Distributions {

private:
    boost::mt19937 rng;

public:

    void set_prng(unsigned int seed) {
        rng = boost::mt19937(seed);
    }

    unsigned int get_random_number() {
        return rng();
    }
 
    boost::mt19937& get_rng() {
        return rng;
    }

    double inv_scaled_chisq_rng(const double a, const double b) {
        return inv_gamma_rng(0.5 * a, 0.5 * a * b);
    }

    double inv_gamma_rng(const double a, const double b) {
        return 1.0 / rgamma(a, 1.0 / b);
    }

    double rgamma(const double a, const double b) {
        boost::random::gamma_distribution<double> myGamma(a, b);
        boost::random::variate_generator<boost::mt19937&, boost::random::gamma_distribution<> > rand_gamma(rng, myGamma);
        double val = rand_gamma();
        return val;
    }

    double beta_rng(const double a, const double b) {
        //std::cout << "@@ beta_rng " << a << ", " << b << std::endl;
        boost::random::beta_distribution<double> mybeta(a, b);
        boost::random::variate_generator<boost::mt19937&, boost::random::beta_distribution<> > rand_beta(rng, mybeta);
        double val = rand_beta();
        //std::cout << "@@ beta_rng val = " << val << std::endl;
        return val;
    }

    double norm_rng(double mean, double sigma2) {
        //std::cout << "@@ norm_rng on " << mean << ", " << sigma2 << std::endl;
        boost::random::normal_distribution<double> nd(mean, std::sqrt(sigma2));
        boost::random::variate_generator< boost::mt19937&, boost::normal_distribution<> > var_nor(rng, nd);
        return var_nor();
    }

    double unif_rng() {
        boost::random::uniform_real_distribution<double> myU(0,1);
        boost::random::variate_generator<boost::mt19937&, boost::random::uniform_real_distribution<> > real_variate_generator(rng, myU);
        return real_variate_generator();
    }

    // Sample from MVN. Reference - https://people.math.sc.edu/Burkardt/cpp_src/normal_dataset/normal_dataset.html
    //upper triangle Cholesky factor R
    double *r8po_fa(int m, std::vector<double> a) {
        double* b;
        int i;
        int j;
        int k;
        double s;

        b = new double[m * m];

        for (j = 0; j < m; j++) {
            for (i = 0; i < m; i++) {
                b[i + j * m] = a[i + j * m];
            }
        }
        for (j = 0; j < m; j++) {
            for (k = 0; k <= j - 1; k++) {
                for (i = 0; i <= k - 1; i++) {
                    b[k + j * m] = b[k + j * m] - b[i + k * m] * b[i + j * m];
                }
                b[k + j * m] = b[k + j * m] / b[k + k * m];
            }
            s = b[j + j * m];
            for (i = 0; i <= j - 1; i++) {
                s = s - b[i + j * m] * b[i + j * m];
            }
            b[j + j * m] = sqrt(s);
        }
        //zero out the lower triangle.
        for (i = 0; i < m; i++) {
            for (j = 0; j < i; j++) {
                b[i + j * m] = 0.0;
            }
        }

        return b;
    }

    double* r8vec_normal_01_new(int t) {
        double* y;
        y = new double[t];
        for (int i = 0; i < t; i++) {
            y[i] = norm_rng(0, 1);
        }
        return y;
    }

    void mvnorm_rng(int m, int n, std::vector<double> a, std::vector<double> mu, int seed, std::vector<double>& mus) {
        int i;
        int j;
        int k;

        //  upper triangular Cholesky factor R of the variance-covariance matrix.
        double *r;
        r = r8po_fa(m, a);

        //  Y = MxN matrix of the 1D normal distribution  mean 0 and variance 1.
        double *y;
        set_prng(seed);
        y = r8vec_normal_01_new(m * n);

        //  Compute X = MU + R' * Y.
        double *x;
        x = new double[m * n];
        for (j = 0; j < n; j++) {
            for (i = 0; i < m; i++) {
                x[i + j * m] = mu[i];
                for (k = 0; k < m; k++) {
                    x[i + j * m] = x[i + j * m] + r[k + i * m] * y[k + j * m];
                }
            }
        }

        delete[] r;
        delete[] y;

        for (j = 0; j < m*n; j++) {
            mus.push_back(x[j]);
        };

        delete[] x;

    }

    // Sample from InvWishart. Reference - https://people.math.sc.edu/Burkardt/cpp_src/wishart/wishart.html
    double* r8mat_copy_new(int m, int n, std::vector<double> a1)

        //  Purpose:
        //
        //    R8MAT_COPY_NEW copies one R8MAT to a "new" R8MAT.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8's, which
        //    may be stored as a vector in column-major order.
        // 
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A1[M*N], the matrix to be copied.
        //
        //    Output, double R8MAT_COPY_NEW[M*N], the copy of A1.
        //
    {
        double* a2;
        int i;
        int j;

        a2 = new double[m * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                a2[i + j * m] = a1[i + j * m];
            }
        }
        return a2;
    }
 
    // The upper Cholesky factor of a symmetric R8MAT
    double* r8mat_cholesky_factor_upper(int n, std::vector<double> a)

        //  Purpose:
        //
        //    R8MAT_CHOLESKY_FACTOR_UPPER: the upper Cholesky factor of a symmetric R8MAT.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    The matrix must be symmetric and positive semidefinite.
        //
        //    For a positive semidefinite symmetric matrix A, the Cholesky factorization
        //    is an upper triangular matrix R such that:
        //
        //      A = R' * R
        //
        //  Parameters:
        //
        //    Input, int N, the number of rows and columns of the matrix A.
        //
        //    Input, double A[N*N], the N by N matrix.
        //    Output, double R8MAT_CHOLESKY_FACTOR[N*N], the N by N upper triangular
        //    Cholesky factor.
    {
        double* c;
        int i;
        int j;
        int k;
        double sum2;


        c = r8mat_copy_new(n, n, a);

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < j; i++)
            {
                c[j + i * n] = 0.0;
            }
            for (i = j; i < n; i++)
            {
                sum2 = c[i + j * n];
                for (k = 0; k < j; k++)
                {
                    sum2 = sum2 - c[k + j * n] * c[k + i * n];
                }
                if (i == j)
                {
                    c[j + i * n] = sqrt(sum2);
                }
                else
                {
                    if (c[j + j * n] != 0.0)
                    {
                        c[j + i * n] = sum2 / c[j + j * n];
                    }
                    else
                    {
                        c[j + i * n] = 0.0;
                    }
                }
            }
        }

        return c;
    }
    
    double* r8ut_inverse(int n, double a[])

        //    R8UT_INVERSE computes the inverse of a R8UT matrix.
        //
        //  Discussion:
        //
        //    The R8UT storage format is used for an M by N upper triangular matrix,
        //    and allocates space even for the zero entries.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N*N], the R8UT matrix.
        //
        //    Output, double R8UT_INVERSE[N*N], the inverse of the upper
        //    triangular matrix.
        //
    {
        double* b;
        int i;
        int j;
        int k;
        b = new double[n * n];

        for (j = n - 1; 0 <= j; j--)
        {
            for (i = n - 1; 0 <= i; i--)
            {
                if (j < i)
                {
                    b[i + j * n] = 0.0;
                }
                else if (i == j)
                {
                    b[i + j * n] = 1.0 / a[i + j * n];
                }
                else if (i < j)
                {
                    b[i + j * n] = 0.0;

                    for (k = i + 1; k <= j; k++)
                    {
                        b[i + j * n] = b[i + j * n] - a[i + k * n] * b[k + j * n];
                    }
                    b[i + j * n] = b[i + j * n] / a[i + i * n];
                }
            }
        }

        return b;
    }

    double* r8mat_mmt_new(int n1, int n2, int n3, double a[], double b[])

        //  Purpose:
        //
        //    R8MAT_MMT_NEW computes C = A * B'.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    For this routine, the result is returned as the function value.
        //
        //  Parameters:
        //
        //    Input, int N1, N2, N3, the order of the matrices.
        //
        //    Input, double A[N1*N2], double B[N3*N2], the matrices to multiply.
        //
        //    Output, double R8MAT_MMT_NEW[N1*N3], the product matrix C = A * B'.
        //
    {
        double* c;
        int i;
        int j;
        int k;

        c = new double[n1 * n3];

        for (i = 0; i < n1; i++)
        {
            for (j = 0; j < n3; j++)
            {
                c[i + j * n1] = 0.0;
                for (k = 0; k < n2; k++)
                {
                    c[i + j * n1] = c[i + j * n1] + a[i + k * n1] * b[j + k * n3];
                }
            }
        }

        return c;
    }
    
    double* wishart_unit_sample_inverse(int m, int df)

        //  Purpose:
        //
        //    WISHART_UNIT_SAMPLE_INVERSE inverts a unit Wishart sample matrix.
        //
        //  Discussion:
        //
        //    This function requires functions from the PDFLIB and RNGLIB libraries.
        //
        //    The "initialize()" function from RNGLIB must be called before using
        //    this function.
        //
        //  Parameters:
        //
        //    Input, int M, the order of the matrix.
        //
        //    Input, int DF, the number of degrees of freedom.
        //    M <= DF.
        //
        //    Output, double WISHART_UNIT_SAMPLE[M*M], the inverse of a
        //    sample matrix from the unit Wishart distribution.
        //
    {
        double* a;
        double* b;
        double* c;
        double df_chi;
        int i;
        int j;

        c = new double[m * m];

        for (i = 0; i < m; i++)
        {
            for (j = 0; j < i; j++)
            {
                c[i + j * m] = 0.0;
            }
            df_chi = (double)(df - i);
            c[i + i * m] = sqrt(r8_chi_sample(df_chi));
            for (j = i + 1; j < m; j++)
            {
                c[i + j * m] = r8_normal_01_sample();
            }
        }

        //for (int i = 0; i < m; i++)
        //{
        //    for (int j = 0; j < m; j++)
        //    {
        //        std::cout << c[j + i * m] << "\n";
        //    };
        //};

        //  Compute B, the inverse of C.
        b = r8ut_inverse(m, c);


        //for (int i = 0; i < m; i++)
        //{
        //    for (int j = 0; j < m; j++)
        //    {
        //        std::cout << b[j + i * m] << "\n";
        //    };
        //};

        //  The inverse of the Wishart sample matrix C'*C is inv(C) * C'.
        a = r8mat_mmt_new(m, m, m, b, b);


        //for (int i = 0; i < m; i++)
        //{
        //    for (int j = 0; j < m; j++)
        //    {
        //        std::cout << a[j + i * m] << "\n";
        //    };
        //};

        //  Free memory.
        delete[] b;
        delete[] c;

        return a;
    }

    double* r8mat_mm_new(int n1, int n2, int n3, double a[], double b[])

        //  Purpose:
        //
        //    R8MAT_MM_NEW multiplies two matrices.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    For this routine, the result is returned as the function value.
        //
        //  Parameters:
        //
        //    Input, int N1, N2, N3, the order of the matrices.
        //
        //    Input, double A[N1*N2], double B[N2*N3], the matrices to multiply.
        //
        //    Output, double R8MAT_MM_NEW[N1*N3], the product matrix C = A * B.
        //
    {
        double* c;
        int i;
        int j;
        int k;

        c = new double[n1 * n3];

        for (i = 0; i < n1; i++)
        {
            for (j = 0; j < n3; j++)
            {
                c[i + j * n1] = 0.0;
                for (k = 0; k < n2; k++)
                {
                    c[i + j * n1] = c[i + j * n1] + a[i + k * n1] * b[k + j * n2];
                }
            }
        }

        return c;
    }
    
    void inv_wishart_rng(int m, int df, std::vector<double> a, std::vector<std::vector<double>>& wis) {
        //double* a;
        double* r;
        double* s;
        double* ua;
        double* uas;

        // Memory allocation to 1D dynamically, size of 1D array will be m*m, created column-wise
        //a = (double*)malloc(m * m * sizeof(double));

        //for (int i = 0; i < m; ++i) {
        //    for (int j = 0; j < m; ++j) {
        //        // mapping 1D array to 2D array column-wise
        //        a[i * m + j] = sigma[j][i];
        //    };
        //};

        //  Get R, the upper triangular Cholesky factor of SIGMA
        r = r8mat_cholesky_factor_upper(m, a);

        //for (int i = 0; i < m; i++)
        //{
        //    for (int j = 0; j < m; j++)
        //    {
        //        std::cout <<  r[j + i * m]<<"\n";
        //    };
        //};


        //  Get S, the inverse of R.
        s = r8ut_inverse(m, r);

        //for (int i = 0; i < m; i++)
        //{
        //    for (int j = 0; j < m; j++)
        //    {
        //        std::cout << s[j + i * m] << "\n";
        //    };
        //};

        //  Get UA, the inverse of a sample from the unit Wishart distribution.
        ua = wishart_unit_sample_inverse(m, df);

        //  Construct the matrix A = S * UA * S'.
        uas = r8mat_mmt_new(m, m, m, ua, s);
        double* inw = r8mat_mm_new(m, m, m, s, uas);

        //  Free memory.
        delete[] r;
        delete[] s;
        delete[] ua;
        delete[] uas;

        for (int i = 0; i < m; i++)
        {
            for (int j = 0; j < m; j++)
            {
                //std::cout << wis.at(j).at(i) << "\n";
                //std::cout << inw[j + i * m] << "\n";
                wis.at(j).at(i) = inw[j + i * m];
            };
        };

    }
};


