#include <iostream>
#include <vector>
#include <string>
#include <numeric>
#include "QAP.h"
#include "Preprocessing.h"

void demo() {
	// Problem instances from QAPLIB
	std::string path = "had12.txt";

	// Initialization of matrices, arrays and other variables
	std::vector<std::vector<int>> F;					// flow matrix
	std::vector<std::vector<int>> D;					// distance matrix
	std::vector<int> initial;						// initial solution (permutation)
	unsigned int size;							// = initial.size();

	// Load the data and process it into F and D matrices
	preprocessing(&path, &F, &D, &size);

	// Initialize the QAP model given F and D matrices
	QAPSolver model(&F, &D, &size);

	// Generate initial solution and solve using the available algorithms
	model.genInitSol("random", &initial, 100);				// generate initial solution with IRG, comment this out if using GA
	model.TS(initial, 50, 7, 10);						// solve using algorithm, in this case Tabu Search
	model.prtSolution();							// print the best permutation, solution and run time

}
