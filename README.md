# pacs-challenge-3

## Description

This program solves the Laplace equation on a two-dimensional grid using both sequential and parallel approaches based on MPI (Message Passing Interface). It computes the numerical solution and compares it with the exact solution, providing the error in the L2 norm. Additionally, it measures and compares the performance of the sequential and parallel methods for different grid sizes.

## Main Features

- **Laplace Equation Solver**: Solves the Laplace equation with Dirichlet boundary conditions on a square grid.
- **Solution Methods**: Implements both a sequential solver and a parallel solver using MPI.
- **Error Calculation**: Compares the numerical solution with the exact solution and calculates the L2 norm of the error.
- **Performance Measurement**: Measures the execution time for both methods (sequential and parallel) on different grid sizes.
- **VTK Output**: Exports the grid solution in a format that can be visualized with ParaView.

## File Structure

- **solver.hpp**: Defines the `Solver` class and its functionalities.
- **solver.cpp**: Implements the `Solver` class, including methods for sequential and parallel solutions.
- **main.cpp**: Contains the `main` function, which manages MPI initialization, runs the solvers, and gathers results.
- **read_Json.hpp**: Header file for reading JSON configuration files.
- **muparser_fun.hpp**: Header file for defining functions using the `muParser` library.
- **Makefile**: Defines the compilation process for the project.
- **README.md**: This file, which describes the project, its features, and usage instructions.

## Requirements

- **C++ Compiler**: Must support C++11 or higher.
- **MPI**: An MPI distribution (e.g., OpenMPI or MPICH) must be installed to run the parallel code.
- [JSON for Modern C++](https://github.com/nlohmann/json) library for JSON file parsing.
- **muParser** library for parsing mathematical expressions and computing derivatives.
- **ParaView**: Software for visualizing VTK output files (optional).

## Compilation and Execution

1. **Compilation**:
   - Compile the project using the provided Makefile:
     ```bash
     make
     ```
2. **Running**:
   - Run the solver:
     ```bash
     mpirun -np <number_of_processes> ./main
     ```
   - Replace `<number_of_processes>` with the desired number of MPI processes.

## VTK Output

The program saves the numerical solution in VTK files within the `output` directory. These files can be opened and visualized using ParaView for visual analysis.

## Input Parameters

The input parameters are provided via a JSON file (`data.json`). The following parameters are supported:

- **grid_dimension**: Specifies the grid dimension for solving the equation.
- **tolerance**: Specifies the tolerance for convergence of the iterative method.
- **max_iterations**: Specifies the maximum number of iterations allowed.
- **force_function**: Specifies the force function in the differential equation.
- **exact_solution**: Specifies the exact solution of the differential equation.
