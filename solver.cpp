#include "solver.h"

solver::solver(int N, int n) {
	//allocate N planets to planet list
	planet_list = new planet;

	//call initialize function
	initialize_planet(N, n);
}

solver::~solver() {
	//delete dynamic allocated variables
	delete[] planet_list;
}

void solver::Euler() {

}

void solver::VV() {

}

void solver::initialize_planet(int N, int n) {
	for (int i = 0; i < N; i++) {
		planet_list[i].set_allocation(n);
	}

	string temp_s;
	double temp_d[9];

	fs.open("planet_initialize.txt");
	//set planet list with file
	for (int i = 0; i < N; i++) {
		//get & set name
		fs >> temp_s;
		planet_list[i].set_planet_name(temp_s);
		
		//get 9 properties ( file must be in order)
		for (int j = 0; j < 9; j++) {
			fs >> temp_d[j];
		}
		//set 9 porperties
		planet_list[i].set_planet_mass(temp_d[0]);
		planet_list[i].set_planet_F_x(temp_d[1], 0);
		planet_list[i].set_planet_F_y(temp_d[2], 0);
		planet_list[i].set_planet_position_x(temp_d[3], 0);
		planet_list[i].set_planet_position_y(temp_d[4], 0);
		planet_list[i].set_planet_v_x(temp_d[5], 0);
		planet_list[i].set_planet_v_y(temp_d[6], 0);
		planet_list[i].set_planet_a_x(temp_d[7], 0); 
		planet_list[i].set_planet_a_y(temp_d[8], 0);
	}
	fs.close();

	

}