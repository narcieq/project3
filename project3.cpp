#include <iostream>
#include <fstream>
#include <iomanip>
#include <cmath>
#include "time_calculation.h"
#include "solver.h"

using namespace std;


int main(int argc, char* argv[])
{
	//get number of planets from command line
	int N = atoi(argv[argc - 1]);

	//get number of time steps from command line
	int n = atoi(argv[argc - 2]);

	//create solver
	solver SOLVER(N);

	//flags
	//-e : use Euler moethod
	//-v : use velocity verlet method
	//-j : add jupiter
	//-u : run unit test

	// check flag
	for (int i = 1; i < argc; i++) {
		//use Euler moethod
		if ((string(argv[i]).find("-") == 0 && string(argv[i]).find("e") != string::npos)) {
			SOLVER.Euler();
		}

		// use velocity verlet method
		else if ((string(argv[i]).find("-") == 0 && string(argv[i]).find("v") != string::npos)) {
			SOLVER.VV();
		}

		//
		else if ((string(argv[i]).find("-") == 0 && string(argv[i]).find("j") != string::npos)) {

		}

		//run unit test
		else if ((string(argv[i]).find("-") == 0 && string(argv[i]).find("u") != string::npos)) {

		}
	}

	


    return 0;
}

