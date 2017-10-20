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
	ofstream fs;
	string filename;
	
	//set file name : 'method'_n_steps_N_planets.txt
	
	if (flag.find("e") == 0) { filename = "Euler_"; }
	else { filename = "VV_"; }
	filename += to_string(n);
	filename += "_steps_";
	filename += to_string(N);
	filename += "_planets.txt";

	cout << filename;

	

	planet* temp_list = SOLVER->get_planet_list();

	fs.open(filename);
	
	int i = 0;
	for (i = 0; i < N; i++) { // print name
		fs << temp_list[i].get_planet_name() << "                                    ";
	}
	fs << endl;

	for (int j = 0; j < n + 1; j++) { // print values of x & y position
		i = 0;
		for (i = 0; i < N; i++) {
			//cout << temp_list[i].get_planet_position_x(j) << "            " << temp_list[i].get_planet_position_y(j) << "            ";
			fs << temp_list[i].get_planet_position_x(j) << "            " << temp_list[i].get_planet_position_y(j) << "            ";
		}
		fs << endl;
		//cout << endl;
		
	}


	fs.close();
}