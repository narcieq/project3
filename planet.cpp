#include "planet.h"

planet::planet() {
	
}

planet::~planet() {
	//delete dynamic allocated variables
	//delete[] F_x;
	//delete[] F_y;
	//delete[] position_x;
	//delete[] position_y;
	//delete[] v_x;
	//delete[] v_y;
	//delete[] a_x;
	//delete[] a_y;
	//delete[] r;
}

//set variable functions n steps
void planet::set_allocation(int n, int number_planets) {
	F_x = new double[n];
	F_y = new double[n];
	position_x = new double[n];
	position_y = new double[n];
	v_x = new double[n];
	v_y = new double[n];
	a_x = new double[n];
	a_y = new double[n];
	r = new double[number_planets - 1];
}

void planet::set_planet_name(string pname) { name = pname; }
void planet::set_planet_mass(double pmass) { mass = pmass; }
void planet::set_planet_F_x(double pF_x, int index) { F_x[index] = pF_x; }
void planet::set_planet_F_y(double pF_y, int index) { F_y[index] = pF_y; }
void planet::set_planet_position_x(double pposiotion_x, int index) { position_x[index] = pposiotion_x; }
void planet::set_planet_position_y(double pposiotion_y, int index) { position_y[index] = pposiotion_y; }
void planet::set_planet_v_x(double pv_x, int index) { v_x[index] = pv_x; }
void planet::set_planet_v_y(double pv_y, int index) { v_y[index] = pv_y; }
void planet::set_planet_a_x(double pa_x, int index) { a_x[index] = pa_x; }
void planet::set_planet_a_y(double pa_y, int index) { a_y[index] = pa_y; }
void planet::set_planet_r(double pr, int index) { r[index] = pr; }

//get varibale functions
string planet::get_planet_name() { return name; }
double planet::get_planet_mass() { return mass; }
double planet::get_planet_F_x(int j) { return F_x[j]; }
double planet::get_planet_F_y(int j) { return F_y[j]; }
double planet::get_planet_position_x(int j) { return position_x[j]; }
double planet::get_planet_position_y(int j) { return position_y[j]; }
double planet::get_planet_v_x(int j) { return v_x[j]; }
double planet::get_planet_v_y(int j) { return v_y[j]; }
double planet::get_planet_a_x(int j) { return a_x[j]; }
double planet::get_planet_a_y(int j) { return a_y[j]; }
double planet::get_planet_r(int index) { return r[index]; }