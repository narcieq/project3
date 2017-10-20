#pragma once
#include "solver.h"
#include <iostream>
using namespace std;

//gravitational constant
const double G = 6.67e-11;
const double pi = acos(-1.0);

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
	double deb = 0;

	for (int i = 0; i < number_planets; i++) {
		for (int j = 1; j < n + 1; j++) {
			//cout << planet_list[i].get_planet_v_x(j - 1) << endl;
			planet_list[i].set_planet_position_x(planet_list[i].get_planet_position_x(j - 1) + h * planet_list[i].get_planet_v_x(j - 1), j);
			deb = planet_list[i].get_planet_position_x(j);
			planet_list[i].set_planet_position_y(planet_list[i].get_planet_position_y(j - 1) + h * planet_list[i].get_planet_v_y(j - 1), j);
			planet_list[i].set_planet_v_x(planet_list[i].get_planet_v_x(j - 1) + h * planet_list[i].get_planet_a_x(j - 1), j);
			planet_list[i].set_planet_v_y(planet_list[i].get_planet_v_y(j - 1) + h * planet_list[i].get_planet_a_y(j - 1), j);
			//distance is calculated and updated in cal_F : x
			planet_list[i].set_planet_F_x(cal_F(j, i, number_planets, 0), j);
			//distance is calculated and updated in cal_F : y
			planet_list[i].set_planet_F_y(cal_F(j, i, number_planets, 1), j);
			planet_list[i].set_planet_a_x(planet_list[i].get_planet_F_x(j) / planet_list[i].get_planet_mass(), j);
			planet_list[i].set_planet_a_y(planet_list[i].get_planet_F_y(j) / planet_list[i].get_planet_mass(), j);
		}
	}
}

void solver::VV(int n) {
	for (int i = 0; i < number_planets; i++) {
		for (int j = 1; j < n + 1; j++) {
			planet_list[i].set_planet_position_x(h * planet_list[i].get_planet_v_x(j - 1) + 0.5 * h * h * planet_list[i].get_planet_a_x(j - 1), j);
			planet_list[i].set_planet_position_y(h * planet_list[i].get_planet_v_y(j - 1) + 0.5 * h * h * planet_list[i].get_planet_a_y(j - 1), j);
			//distance is calculated and updated in cal_F : x
			planet_list[i].set_planet_F_x(cal_F(j, i, number_planets, 0), j);
			//distance is calculated and updated in cal_F : y
			planet_list[i].set_planet_F_y(cal_F(j, i, number_planets, 1), j);
			planet_list[i].set_planet_a_x(planet_list[i].get_planet_F_x(j) / planet_list[i].get_planet_mass(), j);
			planet_list[i].set_planet_a_y(planet_list[i].get_planet_F_y(j) / planet_list[i].get_planet_mass(), j);
			planet_list[i].set_planet_v_x(0.5 * h * (planet_list[i].get_planet_a_x(j - 1) +  planet_list[i].get_planet_a_x(j)), j);
			planet_list[i].set_planet_v_y(0.5 * h * (planet_list[i].get_planet_a_y(j - 1) +  planet_list[i].get_planet_a_y(j)), j);			
		}

		cout << "v : " << planet_list[i].get_planet_v_x(1) << " " << planet_list[i].get_planet_v_y(1) << endl;
		cout << "a : " << planet_list[i].get_planet_a_x(1) << " " << planet_list[i].get_planet_a_y(1) << endl;
		cout << "v : " << planet_list[i].get_planet_v_x(2) << " " << planet_list[i].get_planet_v_y(2) << endl;
		cout << "a : " << planet_list[i].get_planet_a_x(2) << " " << planet_list[i].get_planet_a_y(2) << endl;
	}

}

void solver::initialize_planet(int N, int n) {
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
		
		//calculate 4 properties
		temp = cal_F(0, i, N, 0);
		planet_list[i].set_planet_F_x(temp, 0);
		planet_list[i].set_planet_a_x(temp / planet_list[i].get_planet_mass(), 0);
		temp = cal_F(0, i, N, 1);
		planet_list[i].set_planet_F_y(temp, 0);
		planet_list[i].set_planet_a_y(temp / planet_list[i].get_planet_mass(), 0);
		
	}
	fs.close();

	

}

double solver::cal_F(int step, int index, int number_planets, bool XY) {
	double F = 0;
	double r = 0;
	double x_index, x_i;
	double y_index, y_i;
	for (int i = 0; i < number_planets; i++) {
		if (i != index) {
			x_index = planet_list[index].get_planet_position_x(step);
			x_i = planet_list[i].get_planet_position_x(step);
			y_index = planet_list[index].get_planet_position_y(step);
			y_i = planet_list[i].get_planet_position_y(step);

			r = sqrt(pow((x_index - x_i), 2) + pow((y_index - y_i), 2));
			
			planet_list[i].set_planet_r(r, index);

			if (XY == 0) { //F_x
				F += -(G*planet_list[i].get_planet_mass()*planet_list[index].get_planet_mass()) / pow(r, 3) * (x_index - x_i);
				
			}
			else { //F_y
				F += -(G*planet_list[i].get_planet_mass()*planet_list[index].get_planet_mass()) / pow(r, 3) * (y_index - y_i);
			}
		}
			
		else {
			
			F += 0;
			
		}
	}

	return F;
}

planet* solver::get_planet_list() {
	return planet_list;
}