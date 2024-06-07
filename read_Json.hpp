#include <fstream>
#include <string>
#include <functional>

#include "muparser_fun.hpp"
#include "json.hpp"
using json = nlohmann::json;


void read_Json(const std::string &filename, int &n, double &tol, int &max_iter,
				std::function<double(double, double)> &force, std::function<double(double, double)> &exact) {
	std::ifstream f(filename);
	json data = json::parse(f);

	// Initializing parameters from JSON data
	n = data.value("grid_dimension", 1);
	tol = data.value("tolerance", 1e-3);
    max_iter = data.value("max_iterations", 10);

	// Initializing forcing term
	std::string force_string = data.value("force_function", "");
	MuparserFun force_function(force_string, 2);
	force = force_function;
	
	// Initializing exact solution
	std::string exact_string = data.value("exact_solution", "");
	MuparserFun exact_solution(exact_string, 2);
	exact = exact_solution;

}