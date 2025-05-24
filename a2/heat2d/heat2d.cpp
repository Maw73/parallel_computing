#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>

#include "helpers.hpp"

using namespace std;

int main(int argc, char **argv)
{
    int max_iterations = 1000;
    double epsilon = 1.0e-3;
    bool verify = true, print_config = true;

    //  MPI hint: remember to initialize MPI first 
    int numprocs = 1; // Use MPI process count instead of this

    // Default values for M rows and N columns
    int N = 12;
    int M = 12;

    process_input(argc, argv, N, M, max_iterations, epsilon, verify, print_config);

    if ( print_config )
        std::cout << "Configuration: m: " << M << ", n: " << N << ", max-iterations: " << max_iterations << ", epsilon: " << epsilon << ", processes: " << numprocs << std::endl;

    auto time_1 = chrono::high_resolution_clock::now(); // change to MPI_Wtime() / omp_get_wtime()

    // The main part of the code that needs to use MPI/OpenMP timing routines 
    
    int i, j;
    double diffnorm;
    int iteration_count = 0;

    Mat U(M, N); // MPI: use local sizes with MPI, e.g., recalculate M and N (e.g., M/numprocs + 2)
    Mat W(M, N); // MPI: use local sizes with MPI, e.g., recalculate M and N

    // Init & Boundary (MPI: different)
    for (i = 0; i < M; ++i) {
        for (j = 0; j < N; ++j) {
            W[i][j] = U[i][j] = 0.0;
        }

        W[i][0] = U[i][0] = 0.05; // left side
        W[i][N-1] = U[i][N-1] = 0.1; // right side
    }

    for (j = 0; j < N; ++j) {
        W[0][j] = U[0][j] = 0.02; // top 
        W[M - 1][j] = U[M - 1][j] = 0.2; // bottom 
    }
    // End init

    iteration_count = 0;
    do
    {
        iteration_count++;
        diffnorm = 0.0;

        // MPI Hint: ISend/IRecv the needed rows based on the position - to process above, below or both

        // Compute new values (but not on boundary) 
        // MPI: based on your process you may need to start and stop at the different row index
        for (i = 1; i < M - 1; ++i)
        {
            for (j = 1; j < N - 1; ++j)
            {
                W[i][j] = (U[i][j + 1] + U[i][j - 1] + U[i + 1][j] + U[i - 1][j]) * 0.25;
                diffnorm += (W[i][j] - U[i][j]) * (W[i][j] - U[i][j]);
            }
        }

        // Only transfer the interior points
        // MPI: based on your process you may need to start and stop at the different row index
        for (i = 1; i < M - 1; ++i)
            for (j = 1; j < N - 1; ++j)
                U[i][j] = W[i][j];

        // MPI: make sure that you have the total diffnorm on all processes for exit criteria
        diffnorm = sqrt(diffnorm); // all processes need to know when to stop
        
    } while (epsilon <= diffnorm && iteration_count < max_iterations);
    
    auto time_2 = chrono::high_resolution_clock::now(); // change to MPI_Wtime() / omp_get_wtime()

    // TODO for MPI: collect all local parts of the U matrix, and save it to another "big" matrix
    // that has the same size the whole size: like bigU(M,N) with the whole matrix allocated. 
    // Hint: star with MPI_Gather then extend to MPI_Gatherv

    // Print time measurements 
    cout << "Elapsed time: "; 
    cout << std::fixed << std::setprecision(4) << chrono::duration<double>(time_2 - time_1).count(); // remove for MPI/OpenMP
    // cout << std::fixed << std::setprecision(4) << (time_2 - time_1); // modify accordingly for MPI/OpenMP
    cout << " seconds, iterations: " << iteration_count << endl; 
 
    // Verification (required for MPI)
    if ( verify ) {
        Mat U_sequential(M, N); // init another matrix for the verification

        int iteration_count_seq = 0;
        heat2d_sequential(U_sequential, max_iterations, epsilon, iteration_count_seq); 

        // Here we need both results - from the sequential (U_sequential) and also from the OpenMP/MPI version, then we compare them with the compare(...) function 
        cout << "Verification: " << ( U.compare(U_sequential) && iteration_count == iteration_count_seq ? "OK" : "NOT OK") << std::endl;
    }

    // MPI: do not forget to call MPI_Finalize()
    
    // U.save_to_disk("heat2d.txt"); // not needed

    return 0;
}