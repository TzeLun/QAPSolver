#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <cstdlib> // needed for exit()

void preprocessing(std::string* path, std::vector<std::vector<int>>* F, std::vector<std::vector<int>>* D, unsigned int* size) {
	std::ifstream inFile;
	inFile.open(path->c_str());
	// only proceed if the file successfully opens
	if (inFile.fail()) {
		std::cout << "ERROR::CODE : Fail to open the file, please check if the file exists" << std::endl;
		std::exit(1); // terminate the program
	}
	int val;
	inFile >> val;
	*size = val;
	for (unsigned int i = 0; i < *size; i++) {
		std::vector<int> temp; 
		for (unsigned int j = 0; j < *size; j++) {
			inFile >> val;
			temp.push_back(val);
		}
		F->push_back(temp);
	}
	for (unsigned int i = 0; i < *size; i++) {
		std::vector<int> temp;
		for (unsigned int j = 0; j < *size; j++) {
			inFile >> val;
			temp.push_back(val);
		}
		D->push_back(temp);
	}
	inFile.close();
}
