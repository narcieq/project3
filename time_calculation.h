#pragma once

/*
	class for time calculation
*/

#include <ctime>
#include <chrono>
#include <iostream>
using namespace std;
using namespace std::chrono;

class time_cal {
private:
	high_resolution_clock::time_point start;
	high_resolution_clock::time_point end;
	duration<double> time_span;
public:
	time_cal();
	~time_cal();
	void time_cal_start();
	void time_cal_end();
	duration<double> get_duration();
};