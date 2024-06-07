#include <mpi.h> // MPI library for parallel computing
#include <iostream> // Standard input/output stream library
#include "read_Json.hpp" // Custom function to read data from JSON file
#include "solver.hpp" // Custom class for solving the problem

int main(int argc, char *argv[]) {
    int provided;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided); // Initialize MPI with thread support

    int n; // Number of grid points
    double tol; // Tolerance for convergence
    int max_iter; // Maximum number of iterations
    std::function<double(double, double)> force; // Function representing the forcing term
    std::function<double(double, double)> exact; // Function representing the exact solution

    std::string filename = "data.json"; // JSON file containing problem parameters
    read_Json(filename, n, tol, max_iter, force, exact); // Read data from JSON file

    Solver solver(n, tol, max_iter, force, exact); // Create an instance of the Solver class with the specified parameters

    // Sequential solution
    solver.sequential_solver(); // Solve the problem sequentially without MPI
    
    // Parallel solution
    solver.parallel_solver(); // Solve the problem using parallel computation with MPI

    MPI_Finalize(); // Finalize MPI
    
    return 0; // Exit the program
}
