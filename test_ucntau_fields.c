#include "include/fields_nate.h"
#include <stdio.h>
#include <stdlib.h>

#define MASS_N 1.674927471e-27
#define GRAV 9.80665e0

int main(int argc, char** argv) {
	double x, y, z, fx, fy, fz, u, t, freq;
	x = 0.0;
	y = 0.0;
	z = 0.0;
	t = 0.0;
	freq = 0.0;

	for(int i = 1; i < 100; i++) {
		z = -1.5 + i * (.05/100.0);
//		z = -1.465 + i * (.001/100.0);
		force_(&x, &y, &z, &fx, &fy, &fz, &u, &t, &freq);
		printf("%f, %e\n", z, (u - -2.390245661413933e-26)/(GRAV*MASS_N)*100);
//		printf("%f, %.15e\n", z, u);
	}
	
	return(0);
}