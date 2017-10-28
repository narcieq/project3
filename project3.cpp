#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cmath>
#include "time_calculation.h"
#include "solver.h"


using namespace std;


void output_file(string flag, solver* SOLVER, int N, int n, int final_time);


int main(int argc, char* argv[])
{
	//get number of planets from command line
	int N = atoi(argv[argc - 1]);

	//get number of time steps from command line
	//input of allocation should be n+1
	int n = atoi(argv[argc - 2]);

	int final_time = atoi(argv[argc - 3]);

	//create solver
	//solver SOLVER(N, n, final_time);

	//flags
	//-c : assume the sun is in the center of mass
		//-e : use Euler moethod
		//-v : use velocity verlet method
	//-m : when sun is not in the center of mass
	//-p : perihelion precession of Mercury

	string flag;

	// check flag
	for (int i = 1; i < argc; i++) {
		// when Sun is in the center of mass : -c
		if ((string(argv[i]).find("-") == 0 && string(argv[i]).find("c") != string::npos)){
			//create solver
			solver SOLVER(N, n, final_time);
			for (int k = i; k < argc; k++) {
				//use Euler moethod : -e
				if ((string(argv[k]).find("-") == 0 && string(argv[k]).find("e") != string::npos)) {
					SOLVER.Euler(n);
					flag = 'e';
					output_file(flag, &SOLVER, N, n, final_time);
				}

				// use velocity verlet method : -v
				else if ((string(argv[k]).find("-") == 0 && string(argv[k]).find("v") != string::npos)) {
					SOLVER.VV(n);
					flag = 'v';
					output_file(flag, &SOLVER, N, n, final_time);
				}
			}
		}

		else if ((string(argv[i]).find("-") == 0 && string(argv[i]).find("m") != string::npos)) {
			//create solver
			solver SOLVER(N, n, final_time);
			SOLVER.VV(n, false);//sun is not center of mass
			flag = 'v';
			output_file(flag, &SOLVER, N, n, final_time);
		}
		else if ((string(argv[i]).find("-") == 0 && string(argv[i]).find("p") != string::npos)) {
			/*SOLVER.~solver();
			solver S(N, n + n / 100, final_time);
			bool need_increase = S.check_time_resolution(n + n / 100);
			while (need_increase == false) {
				S.~solver();
				n = n * 10;
				solver S(N, n + n / 100, final_time);
			}
			 S.VV_relative(n + n/100);//sun is not center of mass
			flag = 'p';
			//output_file(flag, &SOLVER, N, n, final_time);

			//double N = S.get_N_ceta();
			double R = S.get_R_ceta();
			cout << "time resolution (step size) : " << 100/n << endl;
			cout << "ceta of perihelion point = " << R << endl;*/


			int N = 2;
			int n =  1e8;
			int n_limit = 1e6; // limitation of steps number for each class
			int limit_repeat;
			int final_time = 11;
			double h_temp;
			double final_limit_time;

			bool need_increase = true;

			double init_x = 0.3075;
			double init_y = 0;
			double init_vx = 0;
			double init_vy = 12.44;
			double init_ax, init_ay;

			planet* temp_list;
			solver temp(true);

			//while (need_increase == true) {
			
			if (n > n_limit) {
			/*	limit_repeat = n / n_limit;
				h_temp = 10.0 / (double)n; // step size 
				cout << h_temp << endl;
				//system("pause");
				final_limit_time = h_temp * n_limit;
	
				for (int i = 0; i < limit_repeat; i++) {
					cout << "   " << i << endl;
					cout << "final_limit_time" << final_limit_time << endl;
					//need have initial value for poistion ,velocity, acceleration.
					temp.init(true, init_x, init_y, init_vx, init_vy, N, n_limit, final_limit_time);
					temp.VV(n_limit);
					temp_list = temp.get_planet_list();
					init_x = temp_list[1].get_planet_position_x(n_limit);
					init_y = temp_list[1].get_planet_position_y(n_limit);
					init_vx = temp_list[1].get_planet_v_x(n_limit);
					init_vy = temp_list[1].get_planet_v_y(n_limit);
					cout << init_x << " last y " << init_y << " last vx " << init_vx << " last vy " << init_vy << endl;
					//system("pause");

				}// 10year end
				cout << init_x << " last y " << init_y << " last vx " << init_vx << " last vy " << init_vy << endl;
				cout << "10 year end "<< endl;
				
				 // 11th year
				int last_period = 3 * n_limit;
				final_limit_time = h_temp * last_period;
				solver temp2(true);
				temp2.init(true, init_x, init_y, init_vx, init_vy, N, last_period, final_limit_time);
				//temp2.VV(last_period);
				cout << "temp2 created" << endl;
				need_increase = temp2.check_time_resolution(last_period);
				cout << "check time resolution end" << endl;
				//if (need_increase == false) {
				cout << "N ceta ('')" << temp2.get_N_ceta() << endl;
			
				*/
				

				//calculate considering relativity
				init_x = 0.3075;
				init_y = 0;
				init_vx = 0;
				init_vy = 12.44;

				limit_repeat = n / n_limit;
				h_temp = 10.0 / (double)n; // step size 
				cout << h_temp << endl;
				//system("pause");
				final_limit_time = h_temp * n_limit;
				
				for (int i = 0; i < limit_repeat; i++) {
					cout << "   " << i << endl;
					//need have initial value for poistion ,velocity, acceleration.
					temp.init(false, init_x, init_y, init_vx, init_vy, N, n_limit, final_limit_time);
					temp.VV_relative(n_limit);
					temp_list = temp.get_planet_list();
					init_x = temp_list[1].get_planet_position_x(n_limit);
					init_y = temp_list[1].get_planet_position_y(n_limit);
					init_vx = temp_list[1].get_planet_v_x(n_limit);
					init_vy = temp_list[1].get_planet_v_y(n_limit);
					cout << init_x << " last y " << init_y << " last vx " << init_vx << " last vy " << init_vy << endl;
					//system("pause");

				}// 10year end
				cout << init_x << " last y " << init_y << " last vx " << init_vx << " last vy " << init_vy << endl;
				cout << "10 year end " << endl;

				// 11th year
				double last_period = 3 * n_limit;
				final_limit_time = h_temp * last_period;
				solver temp2(true);
				temp2.init(false, init_x, init_y, init_vx, init_vy, N, last_period, final_limit_time);
				//temp2.VV(last_period);
				cout << "temp2 created" << endl;
				temp2.VV_relative(last_period);

				double R_tan_perihelion = temp2.return_perihelion(n);
				temp2.cal_R_perihelion(R_tan_perihelion);
				cout << "R_ceta " << temp2.get_R_ceta() << endl;
				
				
				//}
				//else {
				//	n = n * 10;
				//}

			}

		}
	}
    return 0;
}

void output_file(string flag, solver* SOLVER, int N, int n, int final_time) {
	ofstream fs;
	string filename;
	
	//set file name : 'method'_n_steps_N_planets.txt
	
	if (flag.find("e") != string::npos) { filename = "EULER_"; }
	else if(flag.find("v") != string::npos){ filename = "VV_"; }
	else if (flag.find("p") != string::npos) { filename = "Precession_"; }
	filename += to_string(n);
	filename += "_steps_";
	filename += to_string(final_time);
	filename += "yr_";
	filename += to_string(N);
	filename += "_planets.txt";

	cout << filename <<  endl;	

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
				fs << temp_list[i].get_planet_position_x(j) << "            " << temp_list[i].get_planet_position_y(j) << "            ";
			}
			fs << endl;

		}


		fs.close();
}