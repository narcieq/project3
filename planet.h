#pragma once
//planet class header
#include <string>
using namespace std;

class planet {
private:
	string name;
	double mass;

	//will be allocated n elements
	double* F_x;
	double* F_y;
	double* position_x;
	double* position_y;
	double* v_x;
	double* v_y;
	double* a_x;
	double* a_y;
	double* r;

public:
	planet();
	~planet();
	//set functions
	void set_allocation(int n, int number_planets);
	void set_planet_name(string pname);
	void set_planet_mass(double pmass);
	void set_planet_F_x(double pF_x, int index);
	void set_planet_F_y(double pF_y, int index);
	void set_planet_position_x(double pposiotion_x, int index);
	void set_planet_position_y(double pposiotion_y, int index);
	void set_planet_v_x(double pv_y, int index);
	void set_planet_v_y(double pv_y, int index);
	void set_planet_a_x(double pa_y, int index);
	void set_planet_a_y(double pa_y, int index);
	void set_planet_r(double pr, int index);
	//get functions
	string get_planet_name();
	double get_planet_mass();
	double get_planet_F_x(int j);
	double get_planet_F_y(int j);
	double get_planet_position_x(int j);
	double get_planet_position_y(int j);
	double get_planet_v_x(int j);
	double get_planet_v_y(int j);
	double get_planet_a_x(int j);
	double get_planet_a_y(int j);
	double get_planet_r(int index);
};