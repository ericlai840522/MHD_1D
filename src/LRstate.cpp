#include <iostream>
#include <math.h>
#include "struct_vars.h"
#include "LRstate.h"
using namespace std;

void LRstate(struct U_VAR *ULx, struct U_VAR *URx, struct U_VAR *U, int nx){
	// 1st order accuracy
	int im;
	for (int i=0; i<nx; i++){ 
	    im = (i-1+nx)%nx; // im for i minus
		ULx[i] = U[im];
	}

	for(int i=0; i<nx; i++)
		URx[i] = U[i];
	
}
