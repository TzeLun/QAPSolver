/*****************************************************************************************
** Quadratic Assignment Problem Solver							**
** Written by Lok Tze Lun, last updated on 5-1-2022					**
******************************************************************************************/

#include <iostream>
#include <vector>
#include <chrono>
#include <ctime>
#include <cstdlib>
#include <random>
#include <numeric>
#include <cmath>
#include <algorithm>
#pragma once



class QAPSolver {
	public:
		std::vector<std::vector<int>> fMat; //flow matrix, n by n
		std::vector<std::vector<int>> dMat; // distance matrix, n by n
		unsigned int n; // the number of facilities/locations
		QAPSolver(std::vector<std::vector<int>>* flow, std::vector<std::vector<int>>* distance, unsigned int *size) {
			fMat = *flow;
			dMat = *distance;
			n = *size;
		}

		~QAPSolver() {
			std::cout << "\n----QAP SOLVER TERMINATED----\n" << std::endl;
		}

		void genInitSol(const char* method, std::vector<int>* arr, unsigned int iter = 10, const char* mode = "MF2MD") {
			auto start = std::chrono::steady_clock::now();
			if (method == "random") {
				randomGen(&*arr, iter);
			}
			else if (method == "greedy") {
				greedyGen(&*arr, mode);
			}
			else {
				std::cout << "ERROR::CODE -- Unknown Solution Generator Method" << std::endl;
				std::cout << "For generating initial solutions, please select the option of \"random\" or \"greedy\"" << std::endl;
			}
			auto end = std::chrono::steady_clock::now();
			auto diff = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
			runtime_initSol = double(diff.count()) * 0.000001;
		}

		void LocalSearch2Opt(std::vector<int> current, bool firstAccept = false) {
			auto start = std::chrono::steady_clock::now();
			bool breakFlag = false; // to break out from the parent for loop
			bool done = false;
			while (!done) {
				bestSol = current; //make current solution the best
				for (unsigned int i = 0; i < n; i++) {
					for (unsigned int j = i + 1; j < n; j++) {
						//std::cout << "i = " << i << ";" << "j = " << j << std::endl;
						//prtSolution(); // for troubleshooting
						transpose2(current, &i, &j); // current is used to not make any transposition permenant
						if (objectiveFunction(neighborSol) < objectiveFunction(bestSol)) {
							bestSol = neighborSol;
							if (firstAccept == true) {
								breakFlag = true;
								break;
							}
						}
					}
					if (breakFlag) {
						break; // break from parent for loop once first solution is found. Only for First Accept
					}
				}
				if (current == bestSol) {
					done = true; // end the search if no change in solution
				}
				else {
					current = bestSol;
				}
			}
			auto end = std::chrono::steady_clock::now();
			auto diff = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
			runtime = double(diff.count()) * 0.000001;
		}

		void LocalSearch3Opt(std::vector<int> current, bool firstAccept = false) {
			auto start = std::chrono::steady_clock::now();
			bool breakFlag = false; // to break out from the parent for loop
			bool done = false;
			while (!done) {
				bestSol = current; //make current solution the best
				for (unsigned int i = 0; i < n; i++) {
					for (unsigned int j = i + 1; j < n - 1; j++) {
						//std::cout << "i = " << i << ";" << "j = " << j << std::endl;
						//prtSolution(); // for troubleshooting
						transpose3(current, &i, &j); // current is used to not make any transposition permenant
						if (objectiveFunction(neighborSol) < objectiveFunction(bestSol)) {
							bestSol = neighborSol;
							if (firstAccept == true) {
								breakFlag = true;
								break;
							}
						}
					}
					if (breakFlag) {
						break; // break from parent for loop once first solution is found. Only for First Accept
					}
				}
				if (current == bestSol) {
					done = true; // end the search if no change in solution
				}
				else {
					current = bestSol;
				}
			}
			auto end = std::chrono::steady_clock::now();
			auto diff = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
			runtime = double(diff.count()) * 0.000001;
		}

		// Local Search Iterative Improvement with dynamic 2-opt/3-opt
		void LocalSearchDynOpt(std::vector<int> current, bool firstAccept = false, const char* mode = "3-opt") {
			auto start = std::chrono::steady_clock::now();
			bool breakFlag = false; // to break out from the parent for loop
			bool done = false;
			int repeat = 1; // count the number of times no new best solution is found
			while (!done) {
				bestSol = current; //make current solution the best
				if (mode == "2-opt") {
					for (unsigned int i = 0; i < n; i++) {
						for (unsigned int j = i + 1; j < n; j++) {
							transpose2(current, &i, &j); // current is used to not make any transposition permenant
							if (objectiveFunction(neighborSol) < objectiveFunction(bestSol)) {
								bestSol = neighborSol;
								if (firstAccept == true) {
									breakFlag = true;
									break;
								}
							}
						}
						if (breakFlag) {
							break; // break from parent for loop once first solution is found. Only for First Accept
						}
					}
				}
				else if (mode == "3-opt") {
					for (unsigned int i = 0; i < n; i++) {
						for (unsigned int j = i + 1; j < n - 1; j++) {
							transpose3(current, &i, &j); // current is used to not make any transposition permenant
							if (objectiveFunction(neighborSol) < objectiveFunction(bestSol)) {
								bestSol = neighborSol;
								if (firstAccept == true) {
									breakFlag = true;
									break;
								}
							}
						}
						if (breakFlag) {
							break; // break from parent for loop once first solution is found. Only for First Accept
						}
					}
				}
				// check if current best solution is still equal to the new best solution
				if (current == bestSol) {
					repeat++; // incrementing the repeat counter if the solutions are the same
				}
				else {
					current = bestSol;
					repeat = 1; // reset the counter to encourage new searches
				}

				// To alternate the neighborhood operator to encourage exploration
				if (repeat == 2) {
					if (mode == "2-opt") {
						mode = "3-opt";
					}
					else if (mode == "3-opt") {
						mode = "2-opt";
					}
				}
				else if (repeat == 3) {
					done = true; // if the changing the operator has no effect to the solution, terminate the search
				}
			}
			auto end = std::chrono::steady_clock::now();
			auto diff = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
			runtime = double(diff.count()) * 0.000001;
		}

		void simulatedAnnealing(std::vector<int> current, int restart = 1, double T = 100.0, double fTemp = 5.0, float alpha = 0.85, const char* reducRule = "geometric", int iterMax = 100) {
			auto start = std::chrono::steady_clock::now();
			std::srand(std::time(nullptr));
			if (reducRule != "geometric" && reducRule != "linear" && reducRule != "exponential") {
				std::cout << "ERROR::CODE -- Unknown reduction rule" << std::endl;
				std::cout << "Available rules : geometric, linear, exponential" << std::endl;
				std::exit(1);
			}

			if (reducRule == "geometric") {
				if (alpha >= 1 && alpha <= 0) {
					std::cout << "WARNING::MESSAGE -- Solver unable to converge" << std::endl;
					std::cout << "For geometric reduction rule, please select alpha in the range between 0 < alpha < 1" << std::endl;
					std::exit(1);
				}
			}
			unsigned int i;
			unsigned int j;
			double iTemp = T; // initial temperature
			bestSol = current; //make current solution the best
			while (restart > 0) {
				//std::cout << "restart number = " << restart << std::endl;
				//prtSolution(); // for troubleshooting
				while (iTemp >= fTemp) {
					int iter = 0;
					// find all the neighbors with 2-opt
					std::vector<std::vector<int>> neighborhood;
					for (unsigned int i = 0; i < n; i++) {
						for (unsigned int j = i + 1; j < n; j++) {
							transpose2(bestSol, &i, &j);
							neighborhood.push_back(neighborSol);
						}
					}
					while (iter < iterMax) {
						unsigned int randInd = std::rand() % neighborhood.size(); // retrieve a neighbor randomly
						int deltaCost = objectiveFunction(bestSol) - objectiveFunction(neighborhood[randInd]);
						if (deltaCost > 0) {
							bestSol = neighborhood[randInd]; // accept transposition as permanant
						}
						else {
							float P = std::exp(deltaCost / iTemp); // metropolis acceptance criteria
							// check if P is larger than rand(0.0, 1.0)
							if (P >= (static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX))) {
								bestSol = neighborhood[randInd]; // accept transposition as permanant
							}
						}
						iter++;
					}
					// Update the current temperature
					if (reducRule == "geometric") {
						iTemp *= alpha;
					}
					else if (reducRule == "linear") {
						iTemp -= alpha;
					}
					else if (reducRule == "exponential") {
						iTemp /= std::exp(alpha);
					}
				}
				iTemp = T;
				restart--;
			}
			auto end = std::chrono::steady_clock::now();
			auto diff = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
			runtime = double(diff.count()) * 0.000001;
		}

		// Genetic Algorithm - Randomized and Unbiased GA
		void GA(int nGen, int size, float mutationRate = 0.4, int crossoverMaxSwitch = 1, int nMutation = 1, unsigned int iter = 1) {
			// GA generates its own initial population. Run time includes this initial generation stage
			auto start = std::chrono::steady_clock::now();
			std::srand(std::time(nullptr));

			std::vector<std::vector<int>> population;  // reserve memory to reduce runtime memory allocation
			std::vector<int> fitness;
			std::vector<int> arr(n);

			// generate the population according to the input population size
			for (unsigned int i = 0; i < n; i++) {
				arr[i] = int(i + 1);
			}
			int population_size = size;
			std::random_shuffle(arr.begin(), arr.end());
			while (population_size > 0) {
				int score = 2000000000;
				int count = iter;
				std::vector<int> temp(n);
				while (count > 0) {
					std::random_shuffle(arr.begin(), arr.end());
					int cscore = objectiveFunction(arr);
					if (cscore < score) {
						score = cscore;
						temp = arr;
					}
					count--;
				}
				population.push_back(temp);
				population_size--;
			}

			// calculation of fitness value for each individual
			for (unsigned int i = 0; i < population.size(); i++) {
				fitness.push_back(objectiveFunction(population[i])); // the problem is here!!!!!!
			}

			// Start of GA
			int gen = 1;
			float P = 1.0; // probability of accepting a child even though it isn't very qualified
			int period = nGen / 5; // three periods with different acceptance probability of the offspring
			while (gen <= nGen) {

				// Parent Selection, randomly selected from the population
				int p1 = std::rand() % size; // Select the first parent
				int p2 = std::rand() % size; // Select the second parent

				// Perform Random Crossover, with crossoverMaxSwitch as the total number of random gene exchange desired
				int nSwitch = crossoverMaxSwitch;
				std::vector<std::vector<int>> child = { population[p1], population[p2] }; // first and second child
				while (nSwitch > 0) {
					int geneLoc = std::rand() % n; // out of n genes, choose 1 randomly
					int c1Gene1 = child[0][geneLoc]; // location of gene 1 in first child
					int c2Gene2 = child[1][geneLoc]; // location of gene 2 in second child (same location as gene 1 from first child)
					int loc_c1Gene2 = int(std::distance(child[0].begin(), std::find(child[0].begin(), child[0].end(), c2Gene2)));
					int loc_c2Gene1 = int(std::distance(child[1].begin(), std::find(child[1].begin(), child[1].end(), c1Gene1)));
					child[0][geneLoc] = c2Gene2;
					child[1][geneLoc] = c1Gene1;
					child[0][loc_c1Gene2] = c1Gene1;
					child[1][loc_c2Gene1] = c2Gene2;
					nSwitch--;
				}

				// Mutation in offsprings, trigger probability defined by mutationRate, mutation is just a single/multiple random swaps, defined by nMutation
				int nMut = nMutation;
				if ((static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX)) >= (1.0 - mutationRate)) {
					while (nMut > 0) {
						for (unsigned int k = 0; k < child.size(); k++) {
							int i = std::rand() % n;
							int j;
							while (true) {
								j = std::rand() % n;
								if (i != j) {
									break;
								}
							}
							// swap the genes randomly within an individual/chromosome
							int temp = child[k][i];
							child[k][i] = child[k][j];
							child[k][j] = temp;
						}
						nMut--;
					}
				}

				// Replacement of individuals in the population
				// Storing the fitness value for child one and two
				std::vector<int> fChild = { objectiveFunction(child[0]), objectiveFunction(child[1]) };

				// adjusting the acceptance probability periodically according to the total number of generations
				if (gen >= period) {
					//std::cout << "Before : " << period << " ; ";
					P -= 0.25;  // encourage exploitation as the trial runs closer to the target generation
					period += nGen / 5;  // four periods of different probability in total
					//std::cout << "After : " << period << std::endl;
				}

				for (unsigned int i = 0; i < fChild.size(); i++) {
					std::vector<int>::iterator worstFitness = std::max_element(fitness.begin(), fitness.end());
					// take in child if it is worthy, those unworthy have higher chance of acceptance during exploration stage, lower during exploitation stage
					if (fChild[i] < *worstFitness) {
						int loc = int(std::distance(fitness.begin(), worstFitness));
						population[loc] = child[i];
						fitness[loc] = fChild[i];
					}
					else {
						if ((static_cast<float>(std::rand()) / static_cast<float>(RAND_MAX)) >= (1.0 - P)) {
							int loc = int(std::distance(fitness.begin(), worstFitness));
							population[loc] = child[i];
							fitness[loc] = fChild[i];
						}
					}
				}

				std::vector<int>::iterator bestFitness = std::min_element(fitness.begin(), fitness.end());
				bestSol = population[int(std::distance(fitness.begin(), bestFitness))];
				gen++;
			}
			auto end = std::chrono::steady_clock::now();
			auto diff = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
			runtime = double(diff.count()) * 0.000001;
		}

		// Tabu Search algorithm. Takes in the initial solution, maximum iteration, min tabu tenure and max tabu tenure
		void TS(std::vector<int> current, int maxIter = 10, int tmin = 7, int tmax = 10) {
			auto start = std::chrono::steady_clock::now();
			std::srand(std::time(nullptr));
			// generate the sequence of tenure values
			std::vector<int> tenureArr;
			generateTenureList(&tenureArr, &tmin, &tmax);
			unsigned int tenureValInd = 0; // to determine the tenure value assigned to each new tabu-active moves
			// initialization step
			bestSol = current;
			std::vector<std::vector<int>> candidateList; // list of neighborhood solution that is tabu-inactive
			std::vector<int> moveValue; // objective function value of the moves
			std::vector<std::vector<unsigned int>> moveList; // stores the 2-opt swap move
			int iter = 1;
			std::vector<std::vector<unsigned int>> tabuList; // stores the best move as a tabu-active move
			std::vector<int> tabuTenureList;

			while (iter <= maxIter) {
				tenureCheck(&tabuList, &tabuTenureList, &iter); // check if the tenure for all tabu-active moves have expired, remove them if true
				for (unsigned int i = 0; i < n; i++) {
					for (unsigned int j = i + 1; j < n; j++) {
						std::vector<unsigned int> mv{ i, j };
						// only consider a move if it is not tabu active, else check if aspiration criteria is met
						if (!tabuCheck(tabuList, mv)) {
							transpose2(current, &i, &j);
							candidateList.push_back(neighborSol);
							moveValue.push_back(objectiveFunction(neighborSol));
							moveList.push_back(mv);
						}
						else {
							transpose2(current, &i, &j);
							int neighborVal = objectiveFunction(neighborSol);
							if (neighborVal < objectiveFunction(bestSol)) {
								// revoke the tabu status of a move that could give better solution than current best solution
								tabuRevoke(&tabuList, &tabuTenureList, mv);
								candidateList.push_back(neighborSol);
								moveValue.push_back(neighborVal);
								moveList.push_back(mv);
							}
						}
					}
				}
				// Select the best candidate, it is the moves that yield the lowest objective function value
				std::vector<int>::iterator bestNeighbor = std::min_element(moveValue.begin(), moveValue.end());
				int bestNeighborPos = int(std::distance(moveValue.begin(), bestNeighbor));

				// Compare the best candidate to the current best solution
				if (moveValue[bestNeighborPos] < objectiveFunction(bestSol)) {
					bestSol = candidateList[bestNeighborPos];
				}

				// best candidate is made as the current solution for the next iteration
				current = candidateList[bestNeighborPos];
				// make the move that yield the best candidate a tabu-active
				tabuAdd(&tabuList, &tabuTenureList, moveList[bestNeighborPos], &tenureArr, &tenureValInd, iter);
				// clear the current iterations candidate list, move list and move values
				candidateList.clear();
				moveValue.clear();
				moveList.clear();
				iter++; // increase the iterations
			}
			auto end = std::chrono::steady_clock::now();
			auto diff = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
			runtime = double(diff.count()) * 0.000001;
		}

		// Print the solution
		void prtSolution() {
			unsigned int ind = 0;
			std::cout << "Permutation : ";
			std::cout << "(";
			for (ind; ind < n - 1; ind++) {
				std::cout << bestSol[ind] << ",";
			}
			std::cout << bestSol[ind] << ")\n";
			std::cout << "Solution : " << objectiveFunction(bestSol) << std::endl;
			std::cout << "Run Time : " << runtime + runtime_initSol << "s" << std::endl;
		}

		std::vector<int> getPerm() {
			return bestSol;
		}

		std::vector<int> getPermInitSol() {
			return initSol;
		}

		int getSolution() {
			return objectiveFunction(bestSol);
		}

		int getSolutionInitSol() {
			return objectiveFunction(initSol);
		}

		double getRunTime() {
			return runtime;
		}

		double getInitSolGenRunTime() {
			return runtime_initSol;
		}

	private:
		double runtime = 0.0;
		double runtime_initSol = 0.0; // runtime for the initial solution
		std::vector<int> neighborSol; // the temporary solution from a neighborhood operation
		std::vector<int> bestSol; // the solution that gives the best value to the optimization problem
		std::vector<int> initSol; // save the initial solution generated

		// generate the initial solution randomly, set iter to non-zero to search iteratively for a better initial solution
		// *arr is the pointer to the array the initial solution is write into. iter is set to 10 by default.
		void randomGen(std::vector<int>* arr, unsigned int iter = 10) {
			std::srand(std::time(nullptr));
			if (!initSol.empty()) {
				initSol.clear();
			}
			// create the array from 1 to n
			for (unsigned int i = 0; i < n; i++) {
				initSol.push_back(int(i + 1));
			}
			std::random_shuffle(initSol.begin(), initSol.end());
			int score = 2000000000;
			std::vector<int> temp = initSol;
			while (iter != 0) {
				std::random_shuffle(initSol.begin(), initSol.end());
				int cscore = objectiveFunction(initSol);
				if (cscore < score) {
					score = cscore;
					temp = initSol;
				}
				iter--;
			}
			*arr = temp;
			initSol = temp;
		}

		// Greedy Search Algorithm to generate initial solution
		// Default search mode is Minimum Flow to Maximum Distance, "MF2MD"
		void greedyGen(std::vector<int>* arr, const char* mode = "MF2MD") {
			if (mode == "MF2MD") {
				std::vector<int> flowSum;
				std::vector<int> distSum;
				std::vector<int> temp(n, 0); // create a vector of n zeros
				for (unsigned int i = 0; i < fMat.size(); i++) {
					int total = 0;
					for (unsigned int j = 0; j < fMat[i].size(); j++) {
						total += fMat[i][j];
					}
					flowSum.push_back(total);
				}
				for (unsigned int i = 0; i < dMat.size(); i++) {
					int total = 0;
					for (unsigned int j = 0; j < dMat[i].size(); j++) {
						total += dMat[i][j];
					}
					distSum.push_back(total);
				}
				std::vector<int> flowSumRef = flowSum; // mutable reference for flow sum values
				std::vector<int> distSumRef = distSum; // mutable reference for distance sum values
				std::vector<int> fminIndDup; // accomodate duplicates of fmin
				std::vector<int> dmaxIndDup; // accomodate duplicates of dmax
				// Assign the facilities with smallest sum of flows to the location with longest sum of distances
				while (std::find(temp.begin(), temp.end(), 0) != temp.end()) {
					std::vector<int>::iterator f_iter = flowSumRef.begin();
					std::vector<int>::iterator d_iter = distSumRef.begin();
					if (!flowSum.empty()) {
						int fmin = *std::min_element(flowSum.begin(), flowSum.end());
						while ((f_iter = std::find(f_iter, flowSumRef.end(), fmin)) != flowSumRef.end()) {
							fminIndDup.push_back(int(std::distance(flowSumRef.begin(), f_iter++)) + 1);
							flowSum.erase(std::find(flowSum.begin(), flowSum.end(), fmin));
						}
					}
					if (!distSum.empty()) {
						int dmax = *std::max_element(distSum.begin(), distSum.end());
						while ((d_iter = std::find(d_iter, distSumRef.end(), dmax)) != distSumRef.end()) {
							dmaxIndDup.push_back(int(std::distance(distSumRef.begin(), d_iter++)));
							distSum.erase(std::find(distSum.begin(), distSum.end(), dmax));
						}
					}
					// Settle the issue of imbalance number of duplicates between location distance and facility flow
					if (fminIndDup.size() == dmaxIndDup.size()) {
						while (!fminIndDup.empty()) {
							std::vector<int>::iterator i = fminIndDup.begin();
							std::vector<int>::iterator j = dmaxIndDup.begin();
							temp[*j] = *i;
							fminIndDup.erase(i);
							dmaxIndDup.erase(j);
						}
					}
					else if (fminIndDup.size() < dmaxIndDup.size()) {
						while (!fminIndDup.empty()) {
							std::vector<int>::iterator i = fminIndDup.begin();
							std::vector<int>::iterator j = dmaxIndDup.begin();
							temp[*j] = *i;
							fminIndDup.erase(i);
							dmaxIndDup.erase(j);
						}
					}
					else if (fminIndDup.size() > dmaxIndDup.size()) {
						while (!dmaxIndDup.empty()) {
							std::vector<int>::iterator i = fminIndDup.begin();
							std::vector<int>::iterator j = dmaxIndDup.begin();
							temp[*j] = *i;
							fminIndDup.erase(i);
							dmaxIndDup.erase(j);
						}
					}
				}
				*arr = temp;
				initSol = temp;
			}
		}

		// calculate the Objective Function value of the QAP
		int objectiveFunction(std::vector<int> solution) {
			int total = 0;
			for (unsigned int i = 0; i < n; i++) {
				for (unsigned int j = 0; j < n; j++) {
					if (i != j) {
						total += fMat[i][j] * dMat[static_cast<unsigned __int64>(solution[i]) - 1][static_cast<unsigned __int64>(solution[j]) - 1];
					} 
				}
			}
			return total;
		}

		// transposition of 2 elements
		void transpose2(std::vector<int> solution, unsigned int* i, unsigned int* j) {
			neighborSol = solution;
			unsigned int temp = neighborSol[*i];
			neighborSol[*i] = neighborSol[*j];
			neighborSol[*j] = temp;
		}

		// transposition of 3 elements
		void transpose3(std::vector<int> solution, unsigned int* i, unsigned int* j) {
			neighborSol = solution;
			// swap the current element i with the current element j
			unsigned int temp = neighborSol[*i];
			neighborSol[*i] = neighborSol[*j];
			neighborSol[*j] = temp;

			// swap the current element i with the current element j + 1
			temp = neighborSol[*i];
			neighborSol[*i] = neighborSol[static_cast<unsigned __int64>(*j) + 1];
			neighborSol[static_cast<unsigned __int64>(*j) + 1] = temp;
		}

		// generate the initial tenure list in ascending order between tmin and tmax (inclusive)
		void generateTenureList(std::vector<int>* tenure_arr, int* tmin, int* tmax) {
			for (int i = *tmin; i <= *tmax; i++) {
				tenure_arr->push_back(i);
			}
		}

		// check if the move is tabu-active. true if it is tabu active, false otherwise
		bool tabuCheck(std::vector<std::vector<unsigned int>> tList, std::vector<unsigned int> move) {
			for (unsigned int i = 0; i < tList.size(); i++) {
				if (tList[i] == move) {
					return true;
				}
			}
			return false;
		}

		// revoke a tabu status based on a move
		void tabuRevoke(std::vector<std::vector<unsigned int>>* tList, std::vector<int>* tenList, std::vector<unsigned int> move) {
			std::vector<std::vector<unsigned int>> t = *tList;
			std::vector<int> ten = *tenList;
			std::vector<unsigned int> moveMirror{ move[1], move[0] };
			std::vector<unsigned int> pos;
			for (unsigned int i = 0; i < tList->size(); i++) {
				if (t[i] == move || t[i] == moveMirror) {
					pos.push_back(i);
				}
				if (pos.size() == 2) {
					break;
				}
			}
			
			tList->erase(tList->begin() + pos[0]);
			tList->erase(tList->begin() + (pos[1] - 1));
			tenList->erase(tenList->begin() + pos[0]);
			tenList->erase(tenList->begin() + (pos[1] - 1));
		}

		// add a new tabu-active move to the tabu list, assigned a different tenure at each iteration using the systematic dynamic tenure principle
		void tabuAdd(std::vector<std::vector<unsigned int>>* tList, std::vector<int>* tenList, std::vector<unsigned int> move, std::vector<int>* tenure_arr, unsigned int* tenureIndex, int iter) {
			std::vector<unsigned int> moveMirror{ move[1], move[0] }; // create the symmetrical opposite of the swap move
			// make a move tabu by appending it and its symmetric opposite to the tabu list
			tList->push_back(move);
			tList->push_back(moveMirror);
			// add the tabu tenure to the new tabu-active move
			// if the end of the sequence for the range of tabu tenure is reached, reshuffle the sequence randomly and start the index at 0 again.
			if (*tenureIndex == tenure_arr->size()) {
				std::random_shuffle(tenure_arr->begin(), tenure_arr->end());
				*tenureIndex = 0;  // restart the sequence
			}
			std::vector<int> ten_arr = *tenure_arr;
			tenList->push_back(iter + ten_arr[*tenureIndex]);
			tenList->push_back(iter + ten_arr[*tenureIndex]); // the symmetrical opposite move also share the same tenure
			*tenureIndex = *tenureIndex + 1;
		}

		// check the tenure of each tabu-active moves, revoke them if the tenure is expired
		void tenureCheck(std::vector<std::vector<unsigned int>>* tList, std::vector<int>* tenList, int* iter) {
			std::vector<int> ten = *tenList;
			std::vector<unsigned int> pos;
			for (unsigned int i = 0; i < ten.size(); i++) {
				if (ten[i] < *iter) {
					pos.push_back(i); // record all expired tabu-active tenure's index in order
				}
			}
			// revoke tabu status for all tabu-active moves with expired tenure
			for (unsigned int j = 0; j < pos.size(); j++) {
				tList->erase(tList->begin() + (pos[j] - j));
				tenList->erase(tenList->begin() + (pos[j] - j));
			}
		}

};
