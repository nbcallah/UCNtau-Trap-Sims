#include "include/fields.h"
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {
	if(argc != 2) {
		printf("Error! Usage: ./test_ucntau_fields posZ\n");
		return(1);
	}
	double posz = atof(argv[1]);
	double pos[3] = {0.0, 0.0, posz};
	printf("%e\n", fieldstrength(pos));
	
	return(0);
}