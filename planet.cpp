#include "planet.h"

planet::planet() {
	
}

planet::~planet() {
	//delete dynamic allocated variables
}

//set variable functions
void planet::set_allocation(int n) {
	F_x = new double[n];
	F_y = new double[n];
	position_x = new double[n];
	position_y = new double[n];
	v_x = new double[n];
	v_y = new double[n];
	a_x = new double[n];
	a_y = new double[n];
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

//get varibale functions
string planet::get_planet_name() { return name; }
double planet::get_planet_mass() { return mass; }
double planet::get_planet_F_x() { return* F_x; }
double planet::get_planet_F_y() { return* F_y; }
double planet::get_planet_position_x() { return* position_x; }
double planet::get_planet_position_y() { return* position_y; }
double planet::get_planet_v_x() { return* v_x; }
double planet::get_planet_v_y() { return* v_y; }
double planet::get_planet_a_x() { return* a_x; }
double planet::get_planet_a_y() { return* a_y; }