#pragma once
//solver class header
#include "planet.h"
#include <cmath>
#include <fstream>
using namespace std;

class solver {
private:
	planet* planet_list;
	fstream fs;
	int number_planets;
	double h; //step size

public:
	solver(int N, int n, int final_time);//need number of planet & steps when constructed
	~solver();

	void Euler(int n);//solve with Euler method
	void VV(int n);//solve with velocity verlet method (Sun is center of mass)
	void VV(int n, bool is_sun_center_mass);//solve with velocity verlet method (Sun is not center of mass)

	void initialize_planet(int N, int n);//set for initialization
	planet* get_planet_list();

	double cal_F_x(int step, int index, int number_planets);
	double cal_F_y(int step, int index, int number_planets);

	double cal_center_x(int step);
	double cal_center_y(int step);

	void cal_Sun_initial();

	double cal_Angular_Momentum(int index, int step);

	//test stability
	double check_circular(int n);
	bool check_conservative(int n);
};
