#include "time_calculation.h"

void time_cal::time_cal_start() {
	start = high_resolution_clock::now();
	
}

void time_cal::time_cal_end() {
	end = high_resolution_clock::now();
	time_span = end - start;
}

duration<double> time_cal::get_duration() {
	cout << "time consumed : " << time_span.count() << " seconds" << endl;
	return time_span;
}

time_cal::time_cal() {}

time_cal::~time_cal() {}