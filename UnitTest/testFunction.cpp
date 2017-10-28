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
*/

TEST_CASE("check ceta") {
	// compare ceta
	cout << "*** check ceta ***" << endl;

	int N = 2;
	int n = 10000000;
	int n_limit = 1000000;
	int limit_repeat;
	int final_time = 11;
	int h_temp;
	int final_limit_time;
	
	bool need_increase = true;

	double init_x = 0.3075;
	double init_y = 0;
	double init_vx = 0;
	double init_vy = 3.408219e-2 * 365;

	planet* temp_list;
	solver temp(true);

	//while (need_increase == true) {
		limit_repeat = n / n_limit;
		h_temp = 10 / n;
		final_limit_time = h_temp * n_limit;
		if (n > n_limit) {
			for (int i = 0; i < limit_repeat; i++) {
				cout << "   " << i << endl;
				temp.init(init_x, init_y, init_vx, init_vy, N, n_limit, final_limit_time);
				temp.VV(n_limit);
				temp_list = temp.get_planet_list();
				init_x = temp_list[1].get_planet_position_x(n_limit);
				init_y = temp_list[1].get_planet_position_y(n_limit);
				init_vx = temp_list[1].get_planet_v_x(n_limit);
				init_vy = temp_list[1].get_planet_v_y(n_limit);
				//temp.~solver();
				//system("pause");
				//delete[] temp_list;
			}// 100year end
			//temp.~solver();

			// 101th year
			int last_period = 2 * n_limit;
			solver temp2(init_x, init_y, init_vx, init_vy, N, last_period, final_limit_time);
			//temp2.VV(last_period);

			need_increase = temp2.check_time_resolution(last_period);
			//if (need_increase == false) {
				cout << temp2.get_N_ceta() << endl;
				temp2.VV_relative(last_period);
				cout << temp2.get_R_ceta();
			//}
			//else {
			//	n = n * 10;
			//}

		}
	//}

}