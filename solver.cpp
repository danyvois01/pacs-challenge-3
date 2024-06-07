#include "solver.hpp" // Include the header file for the Solver class
#include <mpi.h> // MPI library for parallel computing
#include <omp.h> // OpenMP library for parallel computing
#include <cmath> // Math functions library
#include <iostream> // Input/output stream
#include <fstream> // File stream
#include <chrono> // Chrono library for timing

Solver::Solver(int n, double tol, int max_iter, std::function<double(double, double)> force, std::function<double(double, double)> exact)
    : n(n), tol(tol), max_iter(max_iter), h(1.0 / (n - 1)), force(force), exact(exact) {} // Constructor to initialize Solver object with provided parameters

double Solver::calculate_L2_norm(const std::vector<double>& data) {
    double L2_norm = 0.0; // Initialize L2 norm
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double diff = exact(i * h, j * h) - data[i * n + j]; // Calculate the difference between exact and computed solution
            L2_norm += diff * diff; // Accumulate squared difference
        }
    }
    return std::sqrt(h * L2_norm); // Return the square root of the accumulated squared difference
}

void Solver::exchange_boundary(std::vector<double>& U, int local_n, int rank, int size) {
    std::vector<double> send_buffer(n); // Initialize send buffer
    std::vector<double> recv_buffer(n); // Initialize receive buffer

    if (rank < size - 1) { // If the current process is not the last one
        std::copy(U.end() - 2 * n, U.end() - n, send_buffer.begin()); // Copy boundary points to the send buffer
        MPI_Sendrecv(send_buffer.data(), n, MPI_DOUBLE, rank + 1, 0, // Send and receive boundary points
                     recv_buffer.data(), n, MPI_DOUBLE, rank + 1, 1,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        std::copy(recv_buffer.begin(), recv_buffer.end(), U.end() - n); // Copy received boundary points to the solution vector
    }

    if (rank > 0) { // If the current process is not the first one
        std::copy(U.begin() + n, U.begin() + 2 * n, send_buffer.begin()); // Copy boundary points to the send buffer
        MPI_Sendrecv(send_buffer.data(), n, MPI_DOUBLE, rank - 1, 1, // Send and receive boundary points
                     recv_buffer.data(), n, MPI_DOUBLE, rank - 1, 0,
                     MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        std::copy(recv_buffer.begin(), recv_buffer.end(), U.begin()); // Copy received boundary points to the solution vector
    }
}

void Solver::sequential_solver() {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get the rank of the current process

    if (rank == 0) { // If the current process is the root process
        std::vector<double> U(n * n, 0.0); // Initialize solution vector
        std::vector<double> U_new(n * n, 0.0); // Initialize vector for new solution
        bool converged = false; // Initialize convergence flag
        int iter = 0; // Initialize iteration counter

        // Start the timer
        auto start = std::chrono::high_resolution_clock::now();

        while (!converged && iter < max_iter) { // Continue iteration until convergence or maximum iteration reached
            double error = 0.0; // Initialize error

            // Perform iteration over the grid points
            for (int i = 1; i < n - 1; ++i) {
                for (int j = 1; j < n - 1; ++j) {
                    // Compute the new solution using the Jacobi method
                    U_new[i * n + j] = 0.25 * (U[(i - 1) * n + j] + U[(i + 1) * n + j] +
                                            U[i * n + j - 1] + U[i * n + j + 1] +
                                            h * h * force(i * h, j * h));
                    // Update the maximum error
                    error = std::max(error, std::abs(U_new[i * n + j] - U[i * n + j]));
                }
            }

            // Update the solution vector
            U.swap(U_new);
            // Check for convergence
            converged = (error < tol);
            ++iter; // Increment iteration counter
        }

        // Stop the timer
        auto end = std::chrono::high_resolution_clock::now();
        // Compute the elapsed time
        std::chrono::duration<double> elapsed = end - start;

        // Print the results
        std::cout << "\n=== SEQUENTIAL SOLVER RESULTS ===\n"; // Aggiungere una sezione per risultati sequenziali
        std::cout << "Number of iterations (sequential): " << iter << std::endl;
        std::cout << "L2 norm of the error (sequential): " << calculate_L2_norm(U) << std::endl;
        std::cout << "Total execution time (sequential): " << elapsed.count() << " seconds" << std::endl;
        std::cout << "=================================\n";;
        }
}

void Solver::parallel_solver() {
    int size, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Get the number of processes
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // Get the rank of the current process

    int global_n = n - 2; // Compute the global size of the grid
    int local_n = (global_n % size > rank) ? global_n / size + 1 : global_n / size; // Compute the local size of the grid for each process
    local_n += 2; // Add boundary points

    std::vector<double> local_u(local_n * n, 0.0); // Initialize local solution vector
    std::vector<double> local_u_new(local_n * n, 0.0); // Initialize vector for new local solution
    std::vector<double> local_force(local_n * n, 0.0); // Initialize local force vector

    for (int i = 1; i < local_n - 1; ++i) { // Iterate over local grid points
        for (int j = 1; j < n - 1; ++j) { // Iterate over columns
            local_force[i * n + j] = h * h * force((i + rank * (local_n - 2)) * h, j * h); // Compute the force at each grid point
        }
    }

    bool converged = false; // Initialize convergence flag
    int iter = 0; // Initialize iteration counter

    // Start the timer
    auto start = std::chrono::high_resolution_clock::now();

    while (!converged && iter < max_iter) { // Continue iteration until convergence or maximum iteration reached
        double local_error = 0.0; // Initialize local error

        for (int i = 1; i < local_n - 1; ++i) { // Iterate over local grid rows
            for (int j = 1; j < n - 1; ++j) { // Iterate over columns
                // Compute the new local solution using the Jacobi method
                local_u_new[i * n + j] = 0.25 * (local_u[(i - 1) * n + j] + local_u[(i + 1) * n + j] +
                                                local_u[i * n + j - 1] + local_u[i * n + j + 1] +
                                                local_force[i * n + j]);
                // Update the maximum local error
                local_error = std::max(local_error, std::abs(local_u_new[i * n + j] - local_u[i * n + j]));
            }
        }

        local_u.swap(local_u_new); // Update the local solution vector

        double global_error; // Initialize global error
        MPI_Allreduce(&local_error, &global_error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD); // Reduce local error values to find the global maximum error

        converged = (global_error < tol); // Check for convergence
        ++iter; // Increment iteration counter

        if (size > 1) { // If there are more than one processes
            exchange_boundary(local_u, local_n, rank, size); // Exchange boundary points between neighboring processes
        }
    }

    // Stop the timer
    auto end = std::chrono::high_resolution_clock::now();
    // Compute the elapsed time
    std::chrono::duration<double> elapsed = end - start;

    std::vector<int> send_counts(size, 0), send_start_idx(size, 0); // Initialize vectors for sending data counts and starting indices
    std::vector<int> recv_counts(size, 0), recv_start_idx(size, 0); // Initialize vectors for receiving data counts and starting indices
    if (rank == 0) { // If the current process is the root process
        for (int i = 0; i < size; i++) { // Iterate over all processes
            send_counts[i] = (global_n % size > i) ? global_n / size + 1 : global_n / size; // Compute the number of elements to send to each process
            send_counts[i] += 2; // Add boundary points
            send_counts[i] *= n; // Multiply by the grid dimension
            send_start_idx[i] = (i == 0) ? 0 : send_start_idx[i - 1] + send_counts[i - 1] - 2 * n; // Compute the starting index for sending data to each process
        }
        for (int i = 0; i < size; i++) { // Iterate over all processes
            recv_counts[i] = send_counts[i] - 2 * n; // Compute the number of elements to receive from each process
            recv_start_idx[i] = send_start_idx[i] + n; // Compute the starting index for receiving data from each process
        }
    }


    std::vector<double> data(n * n, 0.0); // Initialize data vector to gather the solution
    MPI_Gatherv(local_u.data() + n, (local_n - 2) * n, MPI_DOUBLE, // Gather the local solutions from all processes
                data.data(), recv_counts.data(), recv_start_idx.data(), MPI_DOUBLE,
                0, MPI_COMM_WORLD);

    if (rank == 0) { // If the current process is the root process
        std::cout << "\n=== PARALLEL SOLVER RESULTS ===\n"; // Aggiungere una sezione per risultati paralleli
        std::cout << "Number of iterations: " << iter << std::endl;
        std::cout << "L2 norm of the error: " << calculate_L2_norm(data) << std::endl;
        std::cout << "Total execution time: " << elapsed.count() << " seconds" << std::endl;
        std::cout << "================================\n";
        export_sol(data, "solution.vtk"); // Export the solution to a file
    }
}

void Solver::export_sol(const std::vector<double>& data, const std::string filename) {
    // opens the file
    std::ofstream vtkFile(filename);

    // check if the file was opened
    if (!vtkFile.is_open()) {
        std::cerr << "Error: could not open file " << filename << std::endl;
        return;
    }

    // Write VTK header
    vtkFile << "# vtk DataFile Version 3.0\n";
    vtkFile << "Scalar Field Data\n";
    vtkFile << "ASCII\n"; // file format

    // Write grid data
    vtkFile << "DATASET STRUCTURED_POINTS\n";                                // format of the dataset
    vtkFile << "DIMENSIONS " << n << " " << n << " " << 1 << "\n";           // number of points in each direction
    vtkFile << "ORIGIN " << 0 << " " << 0 << " 0\n";                         // lower-left corner of the structured grid
    vtkFile << "SPACING" << " " << h << " " << h << " " << 1 << "\n";        // spacing between points in each direction
    vtkFile << "POINT_DATA " << (n) * (n) << "\n";                           // number of points

    // Write scalar field data
    vtkFile << "SCALARS scalars double\n"; // description of the scalar field
    vtkFile << "LOOKUP_TABLE default\n";   // color table

    // Write vector field data
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            vtkFile << data[i * n + j] << "\n";
        }
    }
}