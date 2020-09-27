#include <iostream>
#include <math.h>
#include "struct_vars.h"
#include "C_to_P.h"
#include "cfl.h"
using namespace std;

double cfl(struct U_VAR *U, struct W_VAR *W, int nx, double dx, double gamma){
	
	C_to_P(&U[0], &W[0], nx, gamma);
	double max_ux = 0, a_sq = 0, ca_sq = 0, cfx_sq = 0;
	double temp;
	for(int i=0;i<nx;i++){
		a_sq = gamma*W[i].P/W[i].d;
		ca_sq = ( pow(W[i].Bx, 2) + pow(W[i].By, 2) + pow(W[i].Bz, 2) )/W[i].d;
		// fast magnetosonic wave speed, Eq(A9) in Stone and Gardiner, 2008.
		cfx_sq = 0.5*( a_sq + ca_sq + sqrt( pow(a_sq + ca_sq, 2) - 4*a_sq*pow(W[i].Bx, 2)/W[i].d ) );
		temp = fabs(W[i].ux) + sqrt(cfx_sq);
		max_ux = (max_ux > temp) ? max_ux : temp;
	}

	return dx/max_ux;
}
