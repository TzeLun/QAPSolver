/*****************************************************************************************
** Quadratic Assignment Problem Solver							**
** Written by Lok Tze Lun, last updated on 5-1-2022					**
******************************************************************************************/


#include <iostream>
#include <vector>
#include <string>
#include <numeric>
#include "QAP.h"
#include "Preprocessing.h"
#include "Demo.h" // not needed if not using the demo code

void printInitialSolution(std::vector<int> init);
void printBestSolution(std::vector<int> Sol);
void verifyContents(std::vector<std::vector<int>> F, std::vector<std::vector<int>> D, unsigned int size);

int main() {

	// Problem instances from QAPLIB
	std::string path = "had20.txt";

	// Initialization of matrices, arrays and other variables
	std::vector<std::vector<int>> F; // flow matrix
	std::vector<std::vector<int>> D; // distance matrix
	std::vector<int> initial;
	unsigned int size; // = initial.size();

	// Load the data and process it
	preprocessing(&path, &F, &D, &size);

	// Verify if the contents are loaded correctly, not a QAPSolver method
	verifyContents(F, D, size);

	// Initialize the QAP model
	QAPSolver model(&F, &D, &size);

	// tuning and evaluation bench :
	int instances = 100;								// change this if want to run for more iterations to get better averages
	int optimumVal = 6922;								// the value of the global optimum, change this if using different problem instances
	int count = 1;
	double avg_deviation = 0.0;							// average deviation from the optimum solution
	double avg_solution = 0.0;							// average solution after restarting for n times
	double avg_runtime = 0.0;							// average run time including the initial solution generation
	std::vector<int> bestSol;							// store the best solution out of all the solutions from every restart
	int bestVal = 2000000000;							// store the objective function value of the best solution
	int frequency = 0;								// count how frequent the optimum solution is found


	while (count <= instances) {
		model.genInitSol("random", &initial, 100);				// comment this out if using GA
		model.GA(10000, 200, 0.5, 2, 2, 100);
		//model.simulatedAnnealing(initial, 100, 100.0, 0.01, 0.85, "geometric", 50);
		//model.LocalSearchDynOpt(initial, false, "3-opt");
		//model.TS(initial, 500, 20, 30);

		avg_runtime += (model.getRunTime() + model.getInitSolGenRunTime());
		/*Comment the below section out if the code below is uncommented*/
		avg_solution += model.getSolution();
		avg_deviation += std::abs(model.getSolution() - optimumVal) / float(optimumVal);

		if (model.getSolution() < bestVal) {
			bestVal = model.getSolution(); // make the current solution the best
			bestSol = model.getPerm(); // store the permutation of the best solution
		}

		if (model.getSolution() == optimumVal) {
			frequency++;
		}
		/*Section ends here*/

		/**
		// only uncomment this section if you want to examine the performances of the initial solution generator
		avg_solution += model.getSolutionInitSol();
		avg_deviation += std::abs(model.getSolutionInitSol() - optimumVal) / float(optimumVal);

		if (model.getSolutionInitSol() < bestVal) {
			bestVal = model.getSolutionInitSol(); // make the current solution the best
			bestSol = model.getPermInitSol(); // store the permutation of the best solution
		}

		if (model.getSolutionInitSol() == optimumVal) {
			frequency++;
		}
		/*Section ends here*/

		count++;
	}

	std::cout << "Total Run Time : " << avg_runtime << "s\n";

	avg_runtime /= float(instances);
	avg_solution /= float(instances);
	avg_deviation = avg_deviation / float(instances) * 100.0; // in percentage

	std::cout << "Average Run Time : " << avg_runtime << "s\n";
	std::cout << "Average solution : " << avg_solution << "\n";
	std::cout << "Average Deviation from the Global Optimum : " << avg_deviation << "%\n";
	std::cout << "Frequency of reaching Global Optimum : " << frequency <<  " ; " << frequency / float(instances) * 100 << "%" << std::endl;
	printBestSolution(bestSol);
	std::cout << "Best solution (Value) : " << bestVal << std::endl;

	return 0;
}

void printInitialSolution(std::vector<int> init) {
	std::cout << "initial solution : (";
	for (unsigned int i = 0; i < init.size() - 1; i++) {
		std::cout << init[i] << ", ";
	}
	std::cout << init[init.size() - 1] << ")" << std::endl;
}

void printBestSolution(std::vector<int> Sol) {
	std::cout << "Best solution : (";
	for (unsigned int i = 0; i < Sol.size() - 1; i++) {
		std::cout << Sol[i] << ", ";
	}
	std::cout << Sol[Sol.size() - 1] << ")" << std::endl;
}

void verifyContents(std::vector<std::vector<int>> F, std::vector<std::vector<int>> D, unsigned int size) {
	std::cout << size << "\n\n";
	for (unsigned int i = 0; i < size; i++) {
		for (unsigned int j = 0; j < size; j++) {
			std::cout << F[i][j] << ", ";
		}
		std::cout << std::endl;
	}
	std::cout << "\n\n";
	for (unsigned int i = 0; i < size; i++) {
		for (unsigned int j = 0; j < size; j++) {
			std::cout << D[i][j] << ", ";
		}
		std::cout << std::endl;
	}
	std::cout << "\n\n";
}
