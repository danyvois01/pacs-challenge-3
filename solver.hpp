#ifndef SOLVER_HPP // Header guard to prevent multiple inclusion of the header file
#define SOLVER_HPP

#include <vector> // Standard vector library
#include <functional> // Function objects library

class Solver {
public:
    Solver(int n, double tol, int max_iter, std::function<double(double, double)> force, std::function<double(double, double)> exact);
    // Constructor to initialize Solver object with problem parameters and functions
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
