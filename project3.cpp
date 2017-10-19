#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include "time_calculation.h"
#include "solver.h"


using namespace std;

void output_file(string flag, solver* SOLVER, int N, int n);


int main(int argc, char* argv[])
{
	//get number of planets from command line
	int N = atoi(argv[argc - 1]);

	//get number of time steps from command line
	int n = atoi(argv[argc - 2]);

	int final_time = atoi(argv[argc - 3]);

	//create solver
	solver SOLVER(N, n, final_time);

	//flags
	//-e : use Euler moethod
	//-v : use velocity verlet method

	string flag;

	// check flag
	for (int i = 1; i < argc; i++) {
		//use Euler moethod
		if ((string(argv[i]).find("-") == 0 && string(argv[i]).find("e") != string::npos)) {
			SOLVER.Euler(n);
			flag = 'e';
			output_file(flag, &SOLVER, N, n);
		}

		// use velocity verlet method
		else if ((string(argv[i]).find("-") == 0 && string(argv[i]).find("v") != string::npos)) {
			SOLVER.VV(n);
			flag = 'v';
			output_file(flag, &SOLVER, N, n);
		}
	}

	


    return 0;
}

void output_file(string flag, solver* SOLVER, int N, int n) {
	fstream fs;
	string filename;

	//set file name : 'method'_n_steps_N_planets.txt
	if (flag.find("e") == 0) { filename = "Euler_"; }
	else { filename = "VV_"; }
	filename += n;
	filename += "_steps_";
	filename += N;
	filename += "_planets.txt";

	planet* temp_list = SOLVER->get_planet_list();

	fs.open(filename);
	int i = 0;
	while (i != N) { // print name
		fs << temp_list[i].get_planet_name() << "                                    ";
		i++;
	}
	fs << endl;

	for (int j = 0; j < n + 1; j++) { // print values of x & y position
		for (i = 0; i < N; i++) {
			fs << temp_list[i].get_planet_position_x(j) << "   " << temp_list[i].get_planet_position_y(j) << "   ";
		}
		fs << endl;
	}


	fs.close();
}