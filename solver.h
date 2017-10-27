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
	double perihelion_N_ceta;
	double perihelion_R_ceta;

public:
	solver(bool check);
	solver(int N, int n, int final_time);//need number of planet & steps when constructed
	solver(double x, double y, double vx, double vy, int N, int n, int final_time);
	~solver();

	void Euler(int n);//solve with Euler method
	void VV(int n);//solve with velocity verlet method (Sun is center of mass)
	void VV(int n, bool is_sun_center_mass);//solve with velocity verlet method (Sun is not center of mass)
	
	

	void set_2pi();

	void initialize_planet(int N, int n);//set for initialization
	planet* get_planet_list();

	double cal_F_x(int step, int index, int number_planets);
	double cal_F_y(int step, int index, int number_planets);
	double cal_RF_x(int step, int index, int number_planets);
	double cal_RF_y(int step, int index, int number_planets);

	double cal_center_x(int step);
	double cal_center_y(int step);

	double return_perihelion(int n);

	void cal_Sun_initial();

	double cal_Angular_Momentum(int index, int step);

	//Mercury's perihelion pressecion
	void VV_relative(int n);
	bool check_time_resolution(int n);
	double get_N_ceta();
	double get_R_ceta();
	void init(double x, double y, double vx, double vy, int N, int n, int final_time);

	//test stability
	double check_circular(int n);
	bool check_conservative(int n);
};
