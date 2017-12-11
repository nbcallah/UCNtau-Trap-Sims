#include "include/fields_nate.h"
#include <stdio.h>
#include <stdlib.h>

#define MASS_N 1.674927471e-27
#define GRAV 9.80665e0
#define JTONEV 6.2415091e27

int main(int argc, char** argv) {
	double x, y, z, fx, fy, fz, u, t, freq;
	x = 0.05;
	y = 0.1;
	z = 0.0;
	t = 0.0;
	freq = 0.0;

	for(int i = 1; i < 1000; i++) {
		z = -1.5 + i * (.05/1000.0);
//		z = -1.465 + i * (.001/100.0);
		force_(&x, &y, &z, &fx, &fy, &fz, &u, &t, &freq);
//		printf("%f %e\n", 100*(z+1.5), (u - -2.390245661413933e-26)/(GRAV*MASS_N)*100);
		printf("%f %e\n", 100*(z+1.5), (u - -2.390245661413933e-26)*JTONEV);
//		printf("%f, %.15e\n", z, u);
	}
	
	return(0);
}