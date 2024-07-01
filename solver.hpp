#ifndef SOLVER_HPP // Header guard to prevent multiple inclusion of the header file
#define SOLVER_HPP

#include <vector> // Standard vector library
#include <string> // Standard string library
#include <functional> // Function objects library

//@note use doxygen for comemnting classes and functions. It is more standard in C++ community
class Solver {
public:
    Solver(int n, double tol, int max_iter, std::function<double(double, double)> force, std::function<double(double, double)> exact);
    // Constructor to initialize Solver object with problem parameters and functions
    void export_sol(const std::vector<double>& data, const std::string filename); // Method to export data to a file
    void exchange_boundary(std::vector<double>& U, int local_n, int rank, int size);
    void parallel_solver(); // Method to solve the problem using parallel computation
    void sequential_solver(); // Method to solve the problem sequentially
    double calculate_L2_norm(const std::vector<double>& data); // Method to calculate the L2 norm of a vector

private:
    int n; // Number of grid points
    double tol; // Tolerance for convergence
    int max_iter; // Maximum number of iterations
    double h; // Mesh spacing
    std::function<double(double, double)> force; // Function representing the forcing term
    std::function<double(double, double)> exact; // Function representing the exact solution
};

#endif // SOLVER_HPP
