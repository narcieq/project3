#pragma once
#include "solver.h"
#include <iostream>
using namespace std;

//gdefine pi
const double pi = acos(-1.0);

//define G (gravitational constant)
const double G = 4 * pi*pi;

//speed of light in vacuum (AU/year)
double c = 63072;

solver::solver(bool check) {}

void solver::init(bool NR, double x, double y, double vx, double vy, int N, int n, double final_time) {
	// NR : true --> N / false --> R

	number_planets = N;

	//allocate N planets to planet list
	planet_list = new planet[N];

	h = (double)final_time / (double)n;
	
	//system("pause");

	perihelion_N_ceta = 0;
	perihelion_R_ceta = 0;

	for (int i = 0; i < N; i++) {
		planet_list[i].set_allocation(n+1, N);
	}

	string temp_name;
	double temp_mass;
	double temp;

	fs.open("M.txt"); // 
	//set planet list with file
	for (int i = 0; i < N; i++) {
		//get & set name
		fs >> temp_name;
		planet_list[i].set_planet_name(temp_name);
		fs >> temp_mass;
		planet_list[i].set_planet_mass(temp_mass);

	}

	// scale mass : sun = 1
	for (int i = 1; i < N; i++) {
		planet_list[i].set_planet_mass(planet_list[i].get_planet_mass() / planet_list[0].get_planet_mass());
	}
	planet_list[0].set_planet_mass(planet_list[0].get_planet_mass() / planet_list[0].get_planet_mass());


	//initialize Sun
	planet_list[0].set_planet_position_x(0, 0);
	planet_list[0].set_planet_position_y(0, 0);
	planet_list[0].set_planet_v_x(0, 0);
	planet_list[0].set_planet_v_y(0, 0);

	//initialize Mercury
	planet_list[1].set_planet_position_x(x, 0);
	planet_list[1].set_planet_position_y(y, 0);
	planet_list[1].set_planet_v_x(vx, 0);
	planet_list[1].set_planet_v_y(vy, 0);

	//calculate acceleration of Mercury
	if (NR == true) {
		planet_list[1].set_planet_a_x(cal_F_x(0, 1, N) / planet_list[1].get_planet_mass(), 0);
		planet_list[1].set_planet_a_y(cal_F_y(0, 1, N) / planet_list[1].get_planet_mass(), 0);
	}
	else  {
		planet_list[1].set_planet_a_x(cal_RF_x(0, 1, N) / planet_list[1].get_planet_mass(), 0);
		planet_list[1].set_planet_a_y(cal_RF_y(0, 1, N) / planet_list[1].get_planet_mass(), 0);
	}

	cout << "Mercury 0 step" << endl;
	cout << "x = " << planet_list[1].get_planet_position_x(0) << "y = " << planet_list[1].get_planet_position_y(0) << endl;
	cout << "ax = " << planet_list[1].get_planet_a_x(0) << "ay = " << planet_list[1].get_planet_a_y(0) << endl;
	
}

solver::solver(int N, int n, int final_time) {
	number_planets = N;
	
	//allocate N planets to planet list
	planet_list = new planet[N];

	h = (double)final_time / (double)n;

	perihelion_N_ceta = 0;
	perihelion_R_ceta = 0;

	//call initialize function % n+1 is the size of matrix 
	initialize_planet(N, n+1);

}

solver::solver(double x, double y, double vx, double vy, int N, int n, int final_time) {
	number_planets = N;

	//allocate N planets to planet list
	planet_list = new planet[N];

	h = (double)final_time / (double)n;

	perihelion_N_ceta = 0;
	perihelion_R_ceta = 0;

	for (int i = 0; i < N; i++) {
		planet_list[i].set_allocation(n, N);
	}

	string temp_name;
	double temp_mass;
	double temp;

	fs.open("M.txt");
	//set planet list with file
	for (int i = 0; i < N; i++) {
		//get & set name
		fs >> temp_name;
		planet_list[i].set_planet_name(temp_name);
		fs >> temp_mass;
		planet_list[i].set_planet_mass(temp_mass);
	}

	//initialize Sun
	planet_list[0].set_planet_position_x(0, 0);
	planet_list[0].set_planet_position_y(0, 0);
	planet_list[0].set_planet_v_x(0, 0);
	planet_list[0].set_planet_v_y(0, 0);

	//initialize Mercury
	planet_list[1].set_planet_position_x(x, 0);
	planet_list[1].set_planet_position_y(y, 0);
	planet_list[1].set_planet_v_x(vx, 0);
	planet_list[1].set_planet_v_y(vy, 0);


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

		//set_2pi();

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
			planet_list[i].set_planet_position_y(planet_list[i].get_planet_position_y(j - 1) + h * planet_list[i].get_planet_v_y(j - 1) + 0.5 * h * h * planet_list[i].get_planet_a_y(j - 1), j);
		}
		
		planet_list[0].set_planet_position_x(0, j);
		planet_list[0].set_planet_position_y(0, j);
		
		
		//set_2pi();

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
		//cout << "inside VV " << planet_list[1].get_planet_position_x(j) << " y " << planet_list[1].get_planet_position_y(j) << endl;
		//cout << "inside VV " << planet_list[1].get_planet_v_x(j) << " vy " << planet_list[1].get_planet_v_y(j) << endl;
		//cout << "inside VV " << planet_list[1].get_planet_a_x(j) << " ay " << planet_list[1].get_planet_a_y(j) << endl;
		//
		//system("pause");
		
	}
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
	double vx = planet_list[index].get_planet_v_x(step);
	double vy = planet_list[index].get_planet_v_y(step);
	double total_v = sqrt(pow(vx, 2) + pow(vy, 2));
	double r, x, y;
	x = planet_list[1].get_planet_position_x(step);
	y = planet_list[1].get_planet_position_y(step);
	r = sqrt(x*x + y*y);

	double sin_ceta = vy / total_v*x*r - vx / total_v*y / r;

	double L =  r * planet_list[index].get_planet_mass() * total_v * sin_ceta;
	//cout << "L " << L << endl;
	return L;
}

bool solver::check_time_resolution(int n) {
	//calculate with the Newtonian
	for (int j = 1; j < n + 1; j++) {
		planet_list[1].set_planet_position_x(planet_list[1].get_planet_position_x(j - 1) + h * planet_list[1].get_planet_v_x(j - 1) + 0.5 * h * h * planet_list[1].get_planet_a_x(j - 1), j);
		planet_list[1].set_planet_position_y(planet_list[1].get_planet_position_y(j - 1) + h * planet_list[1].get_planet_v_y(j - 1) + 0.5 * h * h * planet_list[1].get_planet_a_y(j - 1), j);
		

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
	}

	double N_tan_perihelion = return_perihelion(n);
	perihelion_N_ceta = atan(N_tan_perihelion) * 180 * 3600 / pi;
	perihelion_N_ceta = abs(perihelion_N_ceta);

	if (perihelion_N_ceta < 43.0e-2) {
		return false;
	}

	else {
		return true;
	}

}

void solver::VV_relative(int n) { // sun = center of mass
								  // calculate center of the mass (sun needs initilal velocity)

	//calculate considering theory of relativity
	for (int j = 1; j < n + 1; j++) {
		planet_list[1].set_planet_position_x(planet_list[1].get_planet_position_x(j - 1) + h * planet_list[1].get_planet_v_x(j - 1) + 0.5 * h * h * planet_list[1].get_planet_a_x(j - 1), j);
		planet_list[1].set_planet_position_y(planet_list[1].get_planet_position_y(j - 1) + h * planet_list[1].get_planet_v_y(j - 1) + 0.5 * h * h * planet_list[1].get_planet_a_y(j - 1), j);

		planet_list[0].set_planet_position_x(0, j);
		planet_list[0].set_planet_position_y(0, j);

		//distance is calculated and updated in cal_F : x
		//get Newtonian velocity
		planet_list[1].set_planet_F_x(cal_F_x(j, 1, number_planets), j);
		planet_list[1].set_planet_a_x(planet_list[1].get_planet_F_x(j) / planet_list[1].get_planet_mass(), j);
		planet_list[1].set_planet_v_x(planet_list[1].get_planet_v_x(j - 1) + 0.5 * h * (planet_list[1].get_planet_a_x(j - 1) + planet_list[1].get_planet_a_x(j)), j);

		//distance is calculated and updated in cal_F : y
		//get Newtonian velocity
		planet_list[1].set_planet_F_y(cal_F_y(j, 1, number_planets), j);
		planet_list[1].set_planet_a_y(planet_list[1].get_planet_F_y(j) / planet_list[1].get_planet_mass(), j);
		planet_list[1].set_planet_v_y(planet_list[1].get_planet_v_y(j - 1) + 0.5 * h * (planet_list[1].get_planet_a_y(j - 1) + planet_list[1].get_planet_a_y(j)), j);


		for (int k = 0; k < number_planets; k++) {
			

			//calculate Relativistic velocity and acceleration of x direction
			planet_list[k].set_planet_F_x(cal_RF_x(j, k, number_planets), j);
			planet_list[k].set_planet_a_x(planet_list[k].get_planet_F_x(j) / planet_list[k].get_planet_mass(), j);
			planet_list[k].set_planet_v_x(planet_list[k].get_planet_v_x(j - 1) + 0.5 * h * (planet_list[k].get_planet_a_x(j - 1) + planet_list[k].get_planet_a_x(j)), j);

			
			//calculate Relativistic velocity and acceleration of y direction
			planet_list[k].set_planet_F_y(cal_RF_y(j, k, number_planets), j);
			planet_list[k].set_planet_a_y(planet_list[k].get_planet_F_y(j) / planet_list[k].get_planet_mass(), j);
			planet_list[k].set_planet_v_y(planet_list[k].get_planet_v_y(j - 1) + 0.5 * h * (planet_list[k].get_planet_a_y(j - 1) + planet_list[k].get_planet_a_y(j)), j);

		}
		cout << planet_list[1].get_planet_position_x(j) << " " << planet_list[1].get_planet_position_y(j) << endl;
		system("pause");
	}
	//double R_tan_perihelion = return_perihelion(n);
	//perihelion_R_ceta = atan(R_tan_perihelion) * 180 * 3600 / pi;
}

void solver::cal_R_perihelion(double R_tan) {
	perihelion_R_ceta = atan(R_tan) * 180 * 3600 / pi;
}

double solver::return_perihelion(int n) {
	// NR : true --> Newtonian , false --> Realtivity
	// n --> n
	double x = planet_list[1].get_planet_position_x(0);
	double y = planet_list[1].get_planet_position_y(0);
	double vx = planet_list[1].get_planet_v_x(0);
	double vy = planet_list[1].get_planet_v_y(0);
	double ax = planet_list[1].get_planet_a_x(0);
	double ay = planet_list[1].get_planet_a_y(0);
	cout << "x = " << x << "  y = " << y << " vx = " << vx << " vy = " << vy << " ax = " << ax << "ay = " << ay << endl;
	double r2 = x*x + y*y;
	double temp;
	int alpha = 1;
	int perihelion_step = 0;
	double total_time = 0;

	while (total_time < 88.0 / 365.0) {
		temp = pow(planet_list[1].get_planet_position_x(alpha),2) + pow(planet_list[1].get_planet_position_y(alpha), 2);
		if (temp < r2) {
			r2 = temp;
			perihelion_step = alpha;
		}
		alpha++;
		total_time += h;
	}

	double tangent = planet_list[1].get_planet_position_y(perihelion_step) / planet_list[1].get_planet_position_x( perihelion_step);
	cout << "periheion_step" << perihelion_step << endl;
	cout << planet_list[1].get_planet_position_x(perihelion_step) << " " << planet_list[1].get_planet_position_y(perihelion_step) << endl;
	return tangent;
}

double solver::cal_RF_x(int step, int index, int number_planets) {
	double F = 0;
	double r = 0;
	double x_index;
	double x_i;
	double y_index;
	double y_i;
	double l;//magnitude of mercuty's angular momentum per unit mass

	for (int i = 0; i < number_planets; i++) {
		if (i != index) {
			l = cal_Angular_Momentum(index, step) / planet_list[index].get_planet_mass();
			//cout << "l in cal_RF" << l << endl;
			//system("pause");
			x_index = planet_list[index].get_planet_position_x(step);
			x_i = planet_list[i].get_planet_position_x(step);
			y_index = planet_list[index].get_planet_position_y(step);
			y_i = planet_list[i].get_planet_position_y(step);

			r = sqrt(pow((x_index - x_i), 2) + pow((y_index - y_i), 2));
			planet_list[index].set_planet_r(r, i);

			F += -(G*planet_list[i].get_planet_mass()*planet_list[index].get_planet_mass()) / pow(r, 3) * (x_index - x_i) *(1 + (3 * l * l) / (r * r * c * c));
		}
		else { F += 0; }
	}
	return F;
}

double solver::cal_RF_y(int step, int index, int number_planets) {
	double F = 0;
	double r = 0;
	double x_index;
	double x_i;
	double y_index;
	double y_i;
	double l;//magnitude of mercuty's angular momentum per unit mass

	for (int i = 0; i < number_planets; i++) {
		if (i != index) {
			l = cal_Angular_Momentum(index, step) / planet_list[i].get_planet_mass();
			x_index = planet_list[index].get_planet_position_x(step);
			x_i = planet_list[i].get_planet_position_x(step);
			y_index = planet_list[index].get_planet_position_y(step);
			y_i = planet_list[i].get_planet_position_y(step);

			r = sqrt(pow((x_index - x_i), 2) + pow((y_index - y_i), 2));
			planet_list[index].set_planet_r(r, i);

			F += -(G*planet_list[i].get_planet_mass()*planet_list[index].get_planet_mass()) / pow(r, 3) * (y_index - y_i) *(1 + (3 * l * l) / (r * r * c * c));
		}
		else { F += 0; }
	}
	return F;
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
	//n --> n+1

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
			//cout << "inside calFX r : " << r << endl;
			planet_list[index].set_planet_r(r, i);

			F += - (G*planet_list[i].get_planet_mass()*planet_list[index].get_planet_mass()) / pow(r, 3) * (x_index - x_i);
			//cout << planet_list[i].get_planet_mass() << " " << planet_list[index].get_planet_mass() << endl;
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



//get functions
planet* solver::get_planet_list() {
	return planet_list;
}
double solver::get_N_ceta() { return perihelion_N_ceta; }
double solver::get_R_ceta() { return perihelion_R_ceta; }



//for unit test & other tests to check the stability of the program
double solver::check_circular(int n) {
	double epsilon = 0;
	double temp;

	double x1, x2, y1, y2, radius1, radius2, radius_max, radius_min;

	

	for (int i = 1; i < number_planets; i++) {
		for (int j = 0; j < n - 1; j++) {
			x1 = planet_list[i].get_planet_position_x(j);
			//x2 = planet_list[i].get_planet_position_x(j);

			y1 = planet_list[i].get_planet_position_y(j);
			//y2 = planet_list[i].get_planet_position_y(j);

			radius1 = sqrt(pow( x1, 2) + pow( y1, 2));
			//radius2 = sqrt(pow( x2, 2) + pow( y2, 2));

			if (j == 0) {
				radius_max = radius1;
				radius_min = radius1;
				radius2 = radius1;
			}
			else {
				if (radius1 >= radius_max) {
					radius_max = radius1;
				}
				else if (radius1 < radius_min) {
					radius_min = radius1;
				}
			}

			temp = abs(radius_max - radius_min);
			if (temp > epsilon) {
				epsilon = temp;
			}
			//system("pause");
		}
	}
	cout << radius_max << endl;
	cout << radius_min << endl;
	cout << log10((radius_max - radius2) / radius2) << endl;


	return epsilon;
}

bool solver::check_conservative(int n) {
	bool is_conserved = true;
	double kineticE = 0;
	double potentialE = 0;
	double totalE_min, totalE_max;
	double v2 = 0;
	double x, y, r;
	double temp = 0;
	double epsilon = 1.0e-13;
	
		for (int j = 0; j < n+1; j++) {
			x = planet_list[1].get_planet_position_x(j);
			y = planet_list[1].get_planet_position_y(j);
			r = sqrt(x*x + y*y);
			v2 = pow(planet_list[1].get_planet_v_x(j), 2) + pow(planet_list[1].get_planet_v_y(j), 2);

			kineticE = 0.5 * planet_list[1].get_planet_mass() * v2;
			potentialE = -G * planet_list[0].get_planet_mass() * planet_list[1].get_planet_mass() / r;
			temp = kineticE + potentialE;

			if (j != 0) {
				if (totalE_max < temp) {
					totalE_max = temp;
				}
				else if (totalE_min > temp) {
					totalE_min = temp;
				}
			}
			else{
				totalE_min = kineticE + potentialE;
				totalE_max = kineticE + potentialE;
			}
		}

		if (totalE_max - totalE_min < epsilon) {
			is_conserved = true;
		}
		else {
			is_conserved = false;
		}
		cout << totalE_max << endl;
		cout << totalE_min << endl;
		cout << totalE_max - totalE_min << endl;

	return is_conserved;
}

void solver::set_2pi() {
	double C = sqrt(2);
	planet_list[1].set_planet_v_y(2 * pi * C, 0);
	planet_list[1].set_planet_v_x(0, 0);
	cout << "initial velocity of Earth = 2 * pi";
	if (C != 1) {
		cout << " * " << C;
	}
	cout << endl << endl;
}