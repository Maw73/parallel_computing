#include <iostream>
#include <iomanip>
#include <chrono>
#include <cmath>
#include <mpi.h>

#include "helpers.hpp"

//Defining two sets of tags to check for different messages
#define TAG_UP 0
#define TAG_DOWN 1

using namespace std;

int main(int argc, char **argv)
{

    int max_iterations = 1000;
    double epsilon = 1.0e-3;
    bool verify = true, print_config = true;

    //  MPI hint: remember to initialize MPI first
    // int numprocs = 1; // Use MPI process count instead of this
    //initializing MPI
    MPI_Init(&argc, &argv);

    //initializing and retrieving values for rank and size
    int rank, numprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //defining two sets of requests to avoid overwriting data
    MPI_Request req_send_up, req_receive_up;
    MPI_Request req_send_down, req_receive_down;

    // Default values for M rows and N columns
    int N = 12;
    int M = 12;

    process_input(argc, argv, N, M, max_iterations, epsilon, verify, print_config);

    if ( print_config )
        std::cout << "Configuration: m: " << M << ", n: " << N << ", max-iterations: " << max_iterations << ", epsilon: " << epsilon << ", processes: " << numprocs << std::endl;

    // auto time_1 = chrono::high_resolution_clock::now(); // change to MPI_Wtime() / omp_get_wtime()
    //waiting time using MPI_Wtime
    const double start = MPI_Wtime();

    // The main part of the code that needs to use MPI/OpenMP timing routines

    int i, j;
    double diffnorm;
    int iteration_count = 0;
    double global_diffnorm;
    //initializing variables for splitting rows
    //computing the number of necessary rows for each process
    int div = M/numprocs;
    //computing the remainders rows
    int remainder = M % numprocs;
    //adding the remainder rows if we are the last process
    int rows_per_proc = div + (rank == (numprocs-1) ? remainder : 0);

    Mat U(rows_per_proc + 2, N); // MPI: use local sizes with MPI, e.g., recalculate M and N (e.g., M/numprocs + 2)
    Mat W(rows_per_proc + 2, N); // MPI: use local sizes with MPI, e.g., recalculate M and N

    // Init & Boundary (MPI: different)
    // for (i = 0; i < M; ++i) {
    //iterate only the number of rows for each processor + 2 ghost rows
    for (i = 0; i < rows_per_proc + 2; ++i) {
        for (j = 0; j < N; ++j) {
            W[i][j] = U[i][j] = 0.0;
        }

        W[i][0] = U[i][0] = 0.05; // left side
        W[i][N-1] = U[i][N-1] = 0.1; // right side
    }

    // for (j = 0; j < N; ++j) {
    //     W[0][j] = U[0][j] = 0.02; // top
    //     W[M - 1][j] = U[M - 1][j] = 0.2; // bottom
    // }
    if(rank == 0)
    //W and U will consider the first row as 1 because we exclude first ghost row
        for (j = 0; j < N; ++j) {
            W[1][j] = U[1][j] = 0.02; // top
        }
    if(rank == numprocs-1)
    //W and U will consider the last row as rows_per_proc because we exclude last ghost row
    //(rows_per_proc+2)-1 rows and a "-1" for excluding the last ghost row
        for (j = 0; j < N; ++j) {
            W[rows_per_proc][j] = U[rows_per_proc][j] = 0.2; // bottom
        }
    // End init

    iteration_count = 0;
    do
    {
        iteration_count++;
        diffnorm = 0.0;
        // MPI Hint: ISend/IRecv the needed rows based on the position - to process above, below or both
        //initializing a start and stop variable which excludes the ghost rows (the rank_first and rank_last will be treated below)
        int start = 1, stop = rows_per_proc+1;
        if(rank == 0){
            start = 2;
            stop = rows_per_proc + 1;
            //sending from the second to last row of P0 to the first row of P1
            MPI_Isend(&U[rows_per_proc][0], N, MPI_DOUBLE, rank+1, TAG_UP, MPI_COMM_WORLD, &req_send_down);
            //receiving on the last row of P0 from the second row of P1
            MPI_Irecv(&U[rows_per_proc + 1][0], N, MPI_DOUBLE, rank+1, TAG_DOWN, MPI_COMM_WORLD, &req_receive_down);
            MPI_Wait(&req_send_down, MPI_STATUS_IGNORE);
            MPI_Wait(&req_receive_down, MPI_STATUS_IGNORE);
        }
        else if(rank==numprocs-1){
            start = 1;
            stop = rows_per_proc;
            //sending from the second row of P_last to the last row of P_last-1
            MPI_Isend(&U[1][0], N, MPI_DOUBLE, rank-1, TAG_DOWN, MPI_COMM_WORLD, &req_send_up);
            //receiving on the first row of P_last from the second to last row of P_last-1
            MPI_Irecv(&U[0][0], N, MPI_DOUBLE, rank-1, TAG_UP, MPI_COMM_WORLD, &req_receive_up);
            MPI_Wait(&req_send_up, MPI_STATUS_IGNORE);
            MPI_Wait(&req_receive_up, MPI_STATUS_IGNORE);
        }
        else{
            //sending from the second row of P_current to the last row of P_previous
            MPI_Isend(&U[1][0], N, MPI_DOUBLE, rank-1, TAG_DOWN, MPI_COMM_WORLD, &req_send_up);
            //receiving on the first row of P_current from the second to last row of P_previous
            MPI_Irecv(&U[0][0], N, MPI_DOUBLE, rank-1, TAG_UP, MPI_COMM_WORLD, &req_receive_up);
            MPI_Wait(&req_send_up, MPI_STATUS_IGNORE);
            MPI_Wait(&req_receive_up, MPI_STATUS_IGNORE);

            //sending from the second to last row of P-current to the first row of P_next
            MPI_Isend(&U[rows_per_proc][0], N, MPI_DOUBLE, rank+1, TAG_UP, MPI_COMM_WORLD, &req_send_down);
            //receiving on the last row of P_current from the second row of P_next
            MPI_Irecv(&U[rows_per_proc + 1][0], N, MPI_DOUBLE, rank+1, TAG_DOWN, MPI_COMM_WORLD, &req_receive_down);
            MPI_Wait(&req_send_down, MPI_STATUS_IGNORE);
            MPI_Wait(&req_receive_down, MPI_STATUS_IGNORE);
        }


        // Compute new values (but not on boundary)
        // MPI: based on your process you may need to start and stop at the different row index
        // for (i = 1; i < M - 1; ++i)
        // {
        //     for (j = 1; j < N - 1; ++j)
        //     {
        //         W[i][j] = (U[i][j + 1] + U[i][j - 1] + U[i + 1][j] + U[i - 1][j]) * 0.25;
        //         diffnorm += (W[i][j] - U[i][j]) * (W[i][j] - U[i][j]);
        //     }
        // }

        for (i = start; i < stop; ++i)
        {
            for (j = 1; j < N - 1; ++j)
            {
                W[i][j] = (U[i][j + 1] + U[i][j - 1] + U[i + 1][j] + U[i - 1][j]) * 0.25;
                diffnorm += (W[i][j] - U[i][j]) * (W[i][j] - U[i][j]);
            }
        }

        // Only transfer the interior points
        // MPI: based on your process you may need to start and stop at the different row index
        for (i = 1; i < rows_per_proc+1; ++i)
            for (j = 1; j < N - 1; ++j)
            //Copying the matrices, so that we can compute the norm
                U[i][j] = W[i][j];

        // MPI: make sure that you have the total diffnorm on all processes for exit criteria
        // diffnorm = sqrt(diffnorm); // all processes need to know when to stop
        MPI_Allreduce(&diffnorm, &global_diffnorm, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        diffnorm = sqrt(global_diffnorm);

    } while (epsilon <= diffnorm && iteration_count < max_iterations);

    // auto time_2 = chrono::high_resolution_clock::now(); // change to MPI_Wtime() / omp_get_wtime()
    const double end = MPI_Wtime();
    //compute the elapsed time
    const double elapsed = end - start;

    // TODO for MPI: collect all local parts of the U matrix, and save it to another "big" matrix
    // that has the same size the whole size: like bigU(M,N) with the whole matrix allocated.
    // Hint: star with MPI_Gather then extend to MPI_Gatherv
    Mat bigU(M, N);
    //number of rows in each process excluding the first and last ghost rows * number of columns
    int sendcount = rows_per_proc * N;
    //initializing receive_counts and receive_displs
    int receive_counts[numprocs];
    int receive_displs[numprocs];

    int offset = 0;
    for (int i = 0; i < numprocs; i++) {
        int no_rows = div + (i == numprocs - 1 ? remainder : 0);
        receive_counts[i] = no_rows * N;
        receive_displs[i] = offset;
        offset += receive_counts[i];
    }

    double gather_time_start = MPI_Wtime();

    MPI_Gatherv(&U[1][0], sendcount, MPI_DOUBLE, &bigU[0][0], receive_counts, receive_displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double gather_time_end = MPI_Wtime();
    double gather_elapsed = gather_time_end - gather_time_start;

    // Print time measurements
    cout << "Elapsed time: ";
    // cout << std::fixed << std::setprecision(4) << chrono::duration<double>(time_2 - time_1).count(); // remove for MPI/OpenMP
    cout << std::fixed << std::setprecision(4) << elapsed; // modify accordingly for MPI/OpenMP
    cout << " seconds, iterations: " << iteration_count << endl;

    // Verification (required for MPI)
    if ( rank==0 & verify ) {
        cout << "Gather time: " << std::fixed << std::setprecision(4) << gather_elapsed << " seconds\n";

        Mat U_sequential(M, N); // init another matrix for the verification

        int iteration_count_seq = 0;
        heat2d_sequential(U_sequential, max_iterations, epsilon, iteration_count_seq);

        // Here we need both results - from the sequential (U_sequential) and also from the OpenMP/MPI version, then we compare them with the compare(...) function
        cout << "Verification: " << ( bigU.compare(U_sequential) && iteration_count == iteration_count_seq ? "OK" : "NOT OK") << std::endl;
    }

    // MPI: do not forget to call MPI_Finalize()
    MPI_Finalize();

    // U.save_to_disk("heat2d.txt"); // not needed

    return 0;
}