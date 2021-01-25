#include <iostream>
#include <iomanip>
#include <vector>
#include <string>

#include "simple_lu_solve.hpp"

using namespace std;

int main()
{

    std::cout << "Example for usage of LU-Solve" << std::endl;

    // new interface build on stl-container to represent vectors
    // and matrices in a flattened scheme : A[i,j] = A_stl[i * N + j]
    // for J -> number of columns

    // initialize
    size_t N = 5;
    double Tol =1e-10;

    std::vector<double> K(N*N);
    std::vector<double> r(N);

    // fill matrix and right hand side (arbitrary, K not singular)
    for (size_t i = 0; i<N; i++)
    {
        K[i*N+i] = 100;
        r[i] = 4.0/3.876;
        for (int j = -1; j<3; j++)
        {
            K[i*N+(i+j)] += 3.0;
        }
    }

    std::cout << "K-Matrix: " << std::endl;
    for (size_t i = 0; i<N; i++)
    {
        for (size_t j = 0; j<N; j++) std::cout << K[i*N+j] << " ";
        std::cout << std::endl;
    }

    std::cout << "r-Vector: " << std::endl;
    for (size_t i = 0; i<N; i++) std::cout << r[i] << " ";
    std::cout << std::endl;
    
    // 1. Perfrom LU decomposition (P vector is used for pivoting)
    std::vector<double> LU(K);
    std::vector<size_t> P(N);
    LUPDecompose(LU.data(), N, Tol, P.data());

    // 2. Check determinant
    double det = LUPDeterminant(LU.data(), P.data(), N);
    std::cout << "detK: " << det << std::endl;
    if (det <= Tol) return -1;

    // 3. Compute inverse
    std::vector<double> Kinv(N*N, 0);
    LUPInvert(LU.data(), P.data(), N, Kinv.data());

    // 4. Test Inverse
    std::vector<double> I(N*N, 0);
    for (size_t i = 0; i<N; i++)
    {
        for (size_t j = 0; j<N; j++)
        {
            I[i*N+j] = 0.0;
            for (size_t k = 0; k<N; k++)
            {
                I[i*N+j] += K[i*N+k] * Kinv[k*N+j];
            }
        }
    }
    std::cout << "K.Kinv = : " << std::setprecision(3) << std::setw(10) << std::endl;
    for (size_t i = 0; i<N; i++)
    {
        for (size_t j = 0; j<N; j++) std::cout << I[i*N+j] << " ";
        std::cout << std::endl;
    }

    // 5. Solve using inverse
    std::vector<double> result_1(N, 0);
    for (size_t i = 0; i<N; i++)
    {
        for (size_t j = 0; j<N; j++) result_1[i] += Kinv[i*N+j] * r[j];
    }
    std::cout << "result = Kinv.r = " << std::endl;
    for (size_t i = 0; i<N; i++) std::cout << result_1[i] << " ";
    std::cout << std::endl;

    // 6. Solve using Gauss
    std::vector<double> result_2(N, 0);
    LUPSolve(LU.data(), P.data(), r.data(), N, result_2.data());
    std::cout << "result = LUSolve = " << std::endl;
    for (size_t i = 0; i<N; i++) std::cout << result_2[i] << " ";
    std::cout << std::endl;



   
}