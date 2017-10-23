#pragma once
#include "solver.h"
#include <iostream>
using namespace std;

//gdefine pi
const double pi = acos(-1.0);
//define G (gravitational constant)
const double G = 4 * pi*pi;

solver::solver(int N, int n, int final_time) {
	number_planets = N;
	
	//allocate N planets to planet list
	planet_list = new planet[N];

	h = (double)final_time / (double)n;

	//call initialize function
	initialize_planet(N, n+1);
}

solver::~solver() {
	//delete dynamic allocated variables
	delete[] planet_list;
}

void solver::Euler(int n) {
	for (int j = 1; j < n + 1; j++) {
		for (int i = 1; i < number_planets; i++) {
			//cout << planet_list[i].get_planet_v_x(j - 1) << endl;
			planet_list[i].set_planet_position_x(planet_list[i].get_planet_position_x(j - 1) + h * planet_list[i].get_planet_v_x(j - 1), j);
			planet_list[i].set_planet_position_y(planet_list[i].get_planet_position_y(j - 1) + h * planet_list[i].get_planet_v_y(j - 1), j);
			planet_list[i].set_planet_v_x(planet_list[i].get_planet_v_x(j - 1) + h * planet_list[i].get_planet_a_x(j - 1), j);
			planet_list[i].set_planet_v_y(planet_list[i].get_planet_v_y(j - 1) + h * planet_list[i].get_planet_a_y(j - 1), j);
		}

		planet_list[0].set_planet_position_x(0, j);
		planet_list[0].set_planet_position_y(0, j);

		for (int k = 0; k < number_planets; k++) {
			//distance is calculated and updated in cal_F : x
			planet_list[k].set_planet_F_x(cal_F_x(j, k, number_planets), j);
			//distance is calculated and updated in cal_F : y
			planet_list[k].set_planet_F_y(cal_F_y(j, k, number_planets), j);
			planet_list[k].set_planet_a_x(planet_list[k].get_planet_F_x(j) / planet_list[k].get_planet_mass(), j);
			planet_list[k].set_planet_a_y(planet_list[k].get_planet_F_y(j) / planet_list[k].get_planet_mass(), j);
		}
		
	}
}

void solver::VV(int n) {
	for (int j = 1; j < n + 1; j++) {
		for (int i = 1; i < number_planets; i++) {
			planet_list[i].set_planet_position_x(planet_list[i].get_planet_position_x(j - 1) + h * planet_list[i].get_planet_v_x(j - 1) + 0.5 * h * h * planet_list[i].get_planet_a_x(j - 1), j);
			//cout << planet_list[i].get_planet_position_x(j-1) << " " << planet_list[i].get_planet_v_x(j - 1) << " " << planet_list[i].get_planet_a_x(j - 1) << endl;
			//system("pause");
			planet_list[i].set_planet_position_y(planet_list[i].get_planet_position_y(j - 1) + h * planet_list[i].get_planet_v_y(j - 1) + 0.5 * h * h * planet_list[i].get_planet_a_y(j - 1), j);


		}

		planet_list[0].set_planet_position_x(0, j);
		planet_list[0].set_planet_position_y(0, j);

		for (int k = 0; k < number_planets; k++) {
			//distance is calculated and updated in cal_F : x
			planet_list[k].set_planet_F_x(cal_F_x(j, k, number_planets), j);
			//distance is calculated and updated in cal_F : y
			planet_list[k].set_planet_F_y(cal_F_y(j, k, number_planets), j);
			planet_list[k].set_planet_a_x(planet_list[k].get_planet_F_x(j) / planet_list[k].get_planet_mass(), j);
			planet_list[k].set_planet_a_y(planet_list[k].get_planet_F_y(j) / planet_list[k].get_planet_mass(), j);
			planet_list[k].set_planet_v_x(planet_list[k].get_planet_v_x(j - 1) + 0.5 * h * (planet_list[k].get_planet_a_x(j - 1) + planet_list[k].get_planet_a_x(j)), j);
			planet_list[k].set_planet_v_y(planet_list[k].get_planet_v_y(j - 1) + 0.5 * h * (planet_list[k].get_planet_a_y(j - 1) + planet_list[k].get_planet_a_y(j)), j);

		}
		
		/*cout << "v : " << planet_list[i].get_planet_v_x(0) << " " << planet_list[i].get_planet_v_y(0) << endl;
		cout << "a : " << planet_list[i].get_planet_a_x(0) << " " << planet_list[i].get_planet_a_y(0) << endl << endl;;
		cout << "v : " << planet_list[i].get_planet_v_x(1) << " " << planet_list[i].get_planet_v_y(1) << endl;
		cout << "a : " << planet_list[i].get_planet_a_x(1) << " " << planet_list[i].get_planet_a_y(1) << endl << endl;
		cout << "v : " << planet_list[i].get_planet_v_x(2) << " " << planet_list[i].get_planet_v_y(2) << endl;
		cout << "a : " << planet_list[i].get_planet_a_x(2) << " " << planet_list[i].get_planet_a_y(2) << endl << endl;*/
	}
	cout << "vv end" << endl;
}

void solver::VV(int n, bool is_sun_center_mass) {
	// calculate center of the mass (sun needs initilal velocity)
	double center_x = 0;
	double center_y = 0;

	center_x = cal_center_x(0);
	center_y = cal_center_y(0);
	cout << "cx " << center_x << endl;
	cout << "cy " << center_y << endl;

	for (int i = 0; i < number_planets; i++) {
		planet_list[i].set_planet_position_x(planet_list[i].get_planet_position_x(0) - center_x, 0);
		planet_list[i].set_planet_position_y(planet_list[i].get_planet_position_y(0) - center_y, 0);
	}

	//give the Sun initial velocity
	cal_Sun_initial();

	for (int j = 1; j < n + 1; j++) {
		for (int i = 0; i < number_planets; i++) {
			planet_list[i].set_planet_position_x(planet_list[i].get_planet_position_x(j - 1) + h * planet_list[i].get_planet_v_x(j - 1) + 0.5 * h * h * planet_list[i].get_planet_a_x(j - 1), j);
			planet_list[i].set_planet_position_y(planet_list[i].get_planet_position_y(j - 1) + h * planet_list[i].get_planet_v_y(j - 1) + 0.5 * h * h * planet_list[i].get_planet_a_y(j - 1), j);
		}

		for (int k = 0; k < number_planets; k++) {
			//distance is calculated and updated in cal_F : x
			planet_list[k].set_planet_F_x(cal_F_x(j, k, number_planets), j);
			//distance is calculated and updated in cal_F : y
			planet_list[k].set_planet_F_y(cal_F_y(j, k, number_planets), j);
			planet_list[k].set_planet_a_x(planet_list[k].get_planet_F_x(j) / planet_list[k].get_planet_mass(), j);
			planet_list[k].set_planet_a_y(planet_list[k].get_planet_F_y(j) / planet_list[k].get_planet_mass(), j);
			planet_list[k].set_planet_v_x(planet_list[k].get_planet_v_x(j - 1) + 0.5 * h * (planet_list[k].get_planet_a_x(j - 1) + planet_list[k].get_planet_a_x(j)), j);
			planet_list[k].set_planet_v_y(planet_list[k].get_planet_v_y(j - 1) + 0.5 * h * (planet_list[k].get_planet_a_y(j - 1) + planet_list[k].get_planet_a_y(j)), j);

		}
	}
	cout << "vv m end" << endl;
}

void solver::cal_Sun_initial() {
	//calculate initial velocity
	double total_L = 0;
	for (int i = 0; i < number_planets; i++) {
		total_L += cal_Angular_Momentum(i, 0);
	}
	double r = sqrt(pow(planet_list[0].get_planet_position_x(0), 2) + pow(planet_list[0].get_planet_position_y(0), 2));
	planet_list[0].set_planet_r(r, 0);

	double temp_v_x = -(total_L * (-planet_list[0].get_planet_position_y(0))) / (pow(planet_list[0].get_planet_r(0),2)*planet_list[0].get_planet_mass());
	double temp_v_y = -(total_L * (planet_list[0].get_planet_position_x(0))) / (pow(planet_list[0].get_planet_r(0), 2)*planet_list[0].get_planet_mass());
	//set initail velocity
	planet_list[0].set_planet_v_x(temp_v_x, 0);
	planet_list[0].set_planet_v_y(temp_v_y, 0);

	//calculate initail acceleration
	double temp_a_x = cal_F_x(0, 0, number_planets);
	double temp_a_y = cal_F_y(0, 0, number_planets);

	//set initial acceleration
	planet_list[0].set_planet_a_x(temp_a_x, 0);
	planet_list[0].set_planet_a_y(temp_a_y, 0);
}

double solver::cal_Angular_Momentum(int index, int step) {
	double total_v = sqrt(pow(planet_list[index].get_planet_v_x(step), 2) + pow(planet_list[index].get_planet_v_y(step), 2));
	double r, x, y;
	x = planet_list[1].get_planet_position_x(step);
	y = planet_list[1].get_planet_position_y(step);
	r = sqrt(x*x + y*y);

	double L =  r * planet_list[index].get_planet_mass() * total_v;
	return L;
}


double solver::cal_center_x(int step) {
	//calculate total mass of every planet
	double total_mass = 0;
	for (int i = 0; i < number_planets; i++) {
		total_mass += planet_list[i].get_planet_mass();
	}

	double center_x = 0;
	double weight = 0;
	for (int i = 0; i < number_planets; i++) {
		weight += planet_list[i].get_planet_mass() * planet_list[i].get_planet_position_x(step);
		cout << weight << endl;
	}
	center_x = weight / total_mass;

	return center_x;
}

double solver::cal_center_y(int step) {
	//calculate total mass of every planet
	double total_mass = 0;
	for (int i = 0; i < number_planets; i++) {
		total_mass += planet_list[i].get_planet_mass();
	}

	double center_y = 0;
	double weight = 0;
	for (int i = 0; i < number_planets; i++) {
		weight += planet_list[i].get_planet_mass() * planet_list[i].get_planet_position_y(step);
	}
	center_y = weight / total_mass;

	return center_y;
}


void solver::initialize_planet(int N, int n) {
	// *** the file should always keep the information of the Sun on the top ***
	// ***else is fine not being in order ***

	for (int i = 0; i < N; i++) {
		planet_list[i].set_allocation(n, N);
	}

	string temp_s;
	double temp_d[5];
	double temp;

	string name;
	cout << "write initialization file name : ";
	cin >> name;
	name += ".txt";

	fs.open(name);
	//set planet list with file
	for (int i = 0; i < N; i++) {
		//get & set name
		fs >> temp_s;
		planet_list[i].set_planet_name(temp_s);
		
		//get 9 properties ( file must be in order)
		for (int j = 0; j < 5; j++) {
			fs >> temp_d[j];
		}
		//set 5 properties
		planet_list[i].set_planet_mass(temp_d[0]);
		planet_list[i].set_planet_position_x(temp_d[1], 0);
		planet_list[i].set_planet_position_y(temp_d[2], 0);
		planet_list[i].set_planet_v_x(temp_d[3]*(double)365, 0);
		planet_list[i].set_planet_v_y(temp_d[4]*(double)365, 0);
		
	}

	// scale mass : sun = 1
	for (int i = 1; i < N; i++) { 
		planet_list[i].set_planet_mass(planet_list[i].get_planet_mass() / planet_list[0].get_planet_mass());
	}
	planet_list[0].set_planet_mass(planet_list[0].get_planet_mass() / planet_list[0].get_planet_mass());

	for (int i = 0; i < N; i++) {
		//calculate 4 properties
		temp = cal_F_x(0, i, N);
		planet_list[i].set_planet_F_x(temp, 0);
		planet_list[i].set_planet_a_x(temp / planet_list[i].get_planet_mass(), 0);
		temp = cal_F_y(0, i, N);
		planet_list[i].set_planet_F_y(temp, 0);
		planet_list[i].set_planet_a_y(temp / planet_list[i].get_planet_mass(), 0);
	}


	fs.close();
}

double solver::cal_F_x(int step, int index, int number_planets) {
	double F = 0;
	double r = 0;
	double x_index;
	double x_i;
	double y_index;
	double y_i;
	for (int i = 0; i < number_planets; i++) {
		if (i != index) {
			x_index = planet_list[index].get_planet_position_x(step);
			x_i = planet_list[i].get_planet_position_x(step);
			y_index = planet_list[index].get_planet_position_y(step);
			y_i = planet_list[i].get_planet_position_y(step);

			r = sqrt(pow((x_index - x_i), 2) + pow((y_index - y_i), 2));
			planet_list[index].set_planet_r(r, i);

			F += - (G*planet_list[i].get_planet_mass()*planet_list[index].get_planet_mass()) / pow(r, 3) * (x_index - x_i);
		}	
		else { F += 0; }
	}
	return F;
}

double solver::cal_F_y(int step, int index, int number_planets) {
	double F = 0;
	double r = 0;
	double x_index;
	double x_i;
	double y_index;
	double y_i;
	for (int i = 0; i < number_planets; i++) {
		if (i != index) {
			x_index = planet_list[index].get_planet_position_x(step);
			x_i = planet_list[i].get_planet_position_x(step);
			//cout << x_index << endl;
			//cout << x_i << endl;
			y_index = planet_list[index].get_planet_position_y(step);
			y_i = planet_list[i].get_planet_position_y(step);

			r = sqrt(pow((x_index - x_i), 2) + pow((y_index - y_i), 2));
			//	planet_list[i].set_planet_r(r, index);
			planet_list[index].set_planet_r(r, i);
			//cout << "r : " << r << endl;
						
			F += -(G*planet_list[i].get_planet_mass()*planet_list[index].get_planet_mass()) / pow(r, 3) * (y_index - y_i);
		}

		else { F += 0; }
	}
	//cout << "Fy : " << planet_list[0].get_planet_F_y(0) << endl;
	//system("pause");

	return F;
}

planet* solver::get_planet_list() {
	return planet_list;
}



//for unit test & other tests to check the stability of the program
double solver::check_circular(int n) {
	double epsilon = 0;
	double temp;

	double x1, x2, x3, y1, y2, y3, radius1, radius2;

	for (int i = 1; i < number_planets; i++) {
		for (int j = 1; j < n - 3; j++) {
			x1 = planet_list[i].get_planet_position_x(j - 1);
			x2 = planet_list[i].get_planet_position_x(j);
			x3 = planet_list[i].get_planet_position_x(j + 1);
			y1 = planet_list[i].get_planet_position_y(j - 1);
			y2 = planet_list[i].get_planet_position_y(j);
			y3 = planet_list[i].get_planet_position_y(j + 1);

			radius1 = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2));
			radius2 = sqrt(pow(x3 - x2, 2) + pow(y3 - y2, 2));

			temp = abs(radius1 - radius2);
			if (temp > epsilon) {
				epsilon = temp;
			}
		}
	}
	return epsilon;
}

bool solver::check_conservative(int n) {
	bool is_conserved = true;
	double kineticE = 0;
	double potentialE = 0;
	double totalE = 0;
	double v2 = 0;
	double x, y, r;
	double temp = 0;
	double epsilon = 1.0e-4;
	
		for (int j = 0; j < n; j++) {
			x = planet_list[1].get_planet_position_x(j);
			y = planet_list[1].get_planet_position_y(j);
			r = sqrt(x*x + y*y);
			v2 = pow(planet_list[1].get_planet_v_x(j), 2) + pow(planet_list[1].get_planet_v_y(j), 2);

			kineticE = 0.5 * planet_list[1].get_planet_mass() * v2;
			potentialE = -G * planet_list[0].get_planet_mass() * planet_list[1].get_planet_mass() / r;
			temp = kineticE + potentialE;
			if (temp != 0) {
				if (abs(temp - totalE) < epsilon) {
					is_conserved = true;
				}
				else {
					is_conserved = false;
					break;
				}
			}
			else{
				totalE = kineticE + potentialE;
			}
		}
	return is_conserved;
}