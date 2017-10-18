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
	void VV();//solve with velocity verlet method

	void initialize_planet(int N, int n);//set for initialization
	planet* get_planet_list();

	double cal_F(int step, int index, int number_planets, bool XY);
};
