//unit test project
#include "catch.hpp"
#include "../project3/solver.h"
#include "../project3/planet.h"
#include "../project3/time_calculation.h"


//unit test
/*TEST_CASE("check_circular") {
	// test whether the orbit is circular
	// by comparing the radius of the orbit for every step
	// with velocity verlet method

	cout << "*** check circular orbit ***" << endl;

	int N = 2;
	int n = 30000; //steps
	int final_time = 3;

	solver SOLVER(N, n, final_time);

	SOLVER.set_2pi();
	SOLVER.VV(n);
	
	double result = SOLVER.check_circular(n);

	//test pass when 1.0e-6
	double epsilon = 1.0e-6;

	planet* t = SOLVER.get_planet_list();
	ofstream fs;
	
	fs.open("escape.txt");
	for (int i = 0; i < n + 1; i++) {
		fs << t[1].get_planet_position_x(i) << " " << t[1].get_planet_position_y(i) << endl;
	}
	fs.close();
	
	REQUIRE((result < epsilon));
}

TEST_CASE("check_conservative") {
	//check 

	cout << "*** check energy conservation ***" << endl;

	int N = 2;
	int n = 10000;
	int final_time = 1;

	solver Solver(N, n, final_time);
	Solver.set_2pi();
	Solver.VV(n);

	double result = Solver.check_conservative(n);

	REQUIRE(result == true);
}

TEST_CASE("cal_Angular_Momentum") {

	cout << "*** check Angular momentum conservation ***" << endl;

	int N = 2;
	int n = 10000;
	int final_time = 1;

	double epsilon = 1.0e-8;
	bool is_conserved = true;

	solver SolveR(N, n, final_time);
	SolveR.set_2pi();
	SolveR.VV(n);
	
	double temp = 0;
	double AngularMomentum = 0;

	AngularMomentum = SolveR.cal_Angular_Momentum(1, 0);

	double max_AM = AngularMomentum;
	double min_AM = AngularMomentum;

	for (int j = 1; j < n+1; j++) {
		temp = SolveR.cal_Angular_Momentum(1, j);
		if (max_AM <= temp) {
			max_AM = temp;
		}
		else if (min_AM > temp) {
			min_AM = temp;
		}
	}

	if (abs(temp - AngularMomentum) < epsilon) {
		is_conserved = true;
		AngularMomentum = temp;
	}
	else {
		is_conserved = false;
	}
	
	REQUIRE(is_conserved == true);
}

TEST_CASE("time") {
	// compare time of VV and Euler
	cout << "*** check time consumed ***" << endl;

	int N = 2;
	int n = 10000;
	int final_time = 1;

	double epsilon = 1.0e-8;
	bool is_conserved = true;

	solver S(N, n, final_time);
	S.set_2pi();
	
	time_cal timer;
	timer.time_cal_start();
	S.VV(n);
	timer.time_cal_end();
	cout << "Velocity Verlet :: ";
	timer.get_duration();
	cout << "Velocity Verlet :: FLOPS = " << 34 * (n+1) << endl;

	time_cal timer2;
	timer2.time_cal_start();
	S.Euler(n);
	timer2.time_cal_end();
	cout << "Euler's method :: ";
	timer2.get_duration();
	cout << "Euler's method :: FLOPS = " << 22 * (n + 1) << endl;
}


TEST_CASE("return_perihelion") {
	int N = 2;
	int n = 100000;
	int final_time = 10;

	solver S(N, n, final_time);
	//file name init

	double epsilon = 1.0e-8;
	double result = S.return_perihelion(n);

	cout << "result " << result << endl;

	REQUIRE(abs(result) < epsilon);
}*/

TEST_CASE("increase_beta") {
	//testing the manipulation of force formula
	//increase r^2 to r^3

	int N = 2;
	int n = 100000;
	int final_time = 1;

	solver S(N, n, final_time);
	//file name init

	double beta = 2.0;
	double F[5] = { 0, 0.5, 0.7, 0.9, 0.999999 };

	//for(int i = 0 ;i < 5; i++) {
		S.VV_increase_beta(n, beta + 1.0);
	//}



	double epsilon = 1.0e-8;

}