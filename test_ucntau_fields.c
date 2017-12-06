#include "include/fields_nate.h"
#include <stdio.h>
#include <stdlib.h>

#define MASS_N 1.674927471e-27
#define GRAV 9.80665e0

int main(int argc, char** argv) {
	double x, y, z, fx, fy, fz, u, t;
	x = 0.0;
	y = 0.0;
	z = 0.0;
	t = 0.0;

	for(int i = 1; i < 100; i++) {
		z = -1.5 + i * (.01/100.0);
		force_(&x, &y, &z, &fx, &fy, &fz, &u, &t);
		printf("%f, %e\n", z, (u + 2.4283243003838247e-26)/(GRAV*MASS_N)*100);
	}
	
	return(0);
}