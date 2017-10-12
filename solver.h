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

public:
	solver(int N);//need number of planet when constructed
	~solver();

	void Euler();//solve with Euler method
	void VV();//solve with velocity verlet method

	void initialize_planet(int N);//set for initialization
};
