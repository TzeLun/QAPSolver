# QAPSolver

A simple and intuitive Quadratic Assignment Problem solver for C++.

## How to use?
Add all the necessary header and cpp files into your Visual Studio project (or any other directory) and copy the the code from the demo.cpp to get started.
Or refer to this code section below:
```C++
#include <iostream>
#include <vector>
#include <string>
#include <numeric>
#include "QAP.h"                                 // Must include!
#include "Preprocessing.h"                       // Not necessary, but highly recommended to use this preprocessing function

void demo() {
	// Problem instances from QAPLIB
	std::string path = "had12.txt";

	// Initialization of matrices, arrays and other variables
	std::vector<std::vector<int>> F;	// flow matrix
	std::vector<std::vector<int>> D;	// distance matrix
	std::vector<int> initial;		// initial solution (permutation)
	unsigned int size;			// = initial.size();

	// Load the data and process it into F and D matrices
	preprocessing(&path, &F, &D, &size);

	// Initialize the QAP model given F and D matrices
	QAPSolver model(&F, &D, &size);

	// Generate initial solution and solve using the available algorithms
	model.genInitSol("random", &initial, 100);	// generate initial solution with IRG, comment this out if using GA
	model.TS(initial, 50, 7, 10);			// solve using algorithm, in this case Tabu Search
	model.prtSolution();				// print the best permutation, solution and run time
}
```

## What algorithm is available?
### To generate initial solution :
#### 1. Iterative Random Generator
- Given a user-defined number of iterations, the algorithm randomly generates a solution that replaces the best solution if it has a smaller cost at every iteration.
```C++
model.genInitSol("random", &initial, 100);  // for random generator, thrid argument is to adjust the number of iterations
```
#### 2. Greedy Generator
- Assigned the facility with the minimum flow to the location with the maximum distance (MFMD)
```C++
model.genInitSol("greedy", &initial);
```
### To solve the QAP :
#### 1. Iterative-Improvement-based Local Search with 2-Opt
```C++
LocalSearch2Opt(initial, false);  // Second argument dictates if the algorithm uses first or best accept as improvement method
```
#### 2. Iterative-Improvement-based Local Search with 3-Opt
```C++
LocalSearch3Opt(initial, false);  // Second argument dictates if the algorithm uses first or best accept as improvement method
```
#### 3. Iterative-Improvement-based Local Search with Dynamic 2/3-Opt
```C++
// Second argument dictates if the algorithm uses first or best accept as improvement method
// Third argument dictates whether the algorithm should start the search with 2-Opt or 3-Opt
LocalSearchDynOpt(initial, false, "3-opt");
```
#### 4. Simulated Annealing with Restarts
Restart SA for number of times using the previous best solution
```C++
// Second Argument : Number of restarts, Third Argument : Initial temperature
// Fourth Argument : Final temperature, Fifth Argument : Annealing Schedule, alpha
// Sixth Argument : reducution rule (currently there are geometric, linear and exponential rule)
// Seventh Argument : number of iteration per temperature level (number of times of randomly picking a neighbour at a temperature) 
simulatedAnnealing(initial, 1, 100.0, 0.01, 0.85, "geometric", 100);
```
#### 5. Random Genetic Algorithm
* Key features:
	* Initial population created randomly using the same iterative random generator
	* Two parents are selected randomly for crossover
	* Each crossover is a single random pair-wise exchange of genes between two parents, however reciprocation of exchange is needed since no repeatability of genes
	* Occurence of mutation is probabilistic, each mutation is a local pair-wise exhange of genes
	* Offspring is added into the population by either:
  		* Offspring has lower cost (better fitness) than the weakest (highest cost) existing chromosomes
  		* Otherwise, offpsring still replace that chromosome at a probability that varies with the progress in the experiment
	* The probability of acceptance starts at 1.0, and decrease by 0.25 at a period of every one-fifth of the total number of generations
<br />
**Note :** When using GA, don't have to generate the initial solution manually!
<br />
```C++
// First argument : number of generations, second argument : population size
// Third argument : mutation rate (probability), Fourth argument : Number of crossover switches per generation
// Fifth argument : number of mutations per generation, Sixth argument : number of iterations (for random generation of initial population)
GA(10000, 50, 0.5, 2, 2, 100);
```
#### 6. Tabu Search
* Key features:
	* Tabu attribute - the move operator (a 2-Opt pair-wise exchange of facilities) and its mirror opposite (ie : swap(2, 3) and swap(3, 2))
	* Aspiration criterion - aspiration by objective function, revoke tabu-active status if tabu moves yield better solution
	* Tabu tenure - random-dynamic tabu tenure, given a range of tenure values, moves are tabu-ed for a tenure assigned sequentially. Once end of sequence is reached, randomly shuffle the tabu tenures to get a new sequence.
<br />
```C++
// Second argument : Number of iterations, number of times the tabu search is run
// Third argument : minimum tabu tenure, Fourth argument : maximum tabu tenure
TS(initial, 10, 7, 10);
```
## Datasets
In general, all dataset can be found in [QAPLIB](https://www.opt.math.tugraz.at/qaplib/inst.html). To get started, you can find the "had12" problem instance in this repository by [Hadley](https://www.opt.math.tugraz.at/qaplib/inst.html#HRW). The format of the problem instances is same as that found on the website.
