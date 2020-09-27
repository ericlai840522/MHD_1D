#include<iostream>
#include<math.h>
#include "struct_vars.h"
#include "LRstate.h"
#include "C_to_P.h"
using namespace std;

void C_to_P(struct U_VAR *U, struct W_VAR *W, int nx, double gamma){
	// Conservative variables to primitive variables
	for (int i=0; i<nx; i++){ 
		W[i].d = U[i].d;
		W[i].ux = U[i].Mx/U[i].d;
		W[i].uy = U[i].My/U[i].d;
		W[i].uz = U[i].Mz/U[i].d;
		W[i].Bx = U[i].Bx;
		W[i].By = U[i].By;
		W[i].Bz = U[i].Bz;
		W[i].P = U[i].E - 0.5*( pow(U[i].Mx, 2) + pow(U[i].My, 2) + pow(U[i].Mz, 2) )/U[i].d - \
					 0.5*( pow(U[i].Bx, 2) + pow(U[i].By, 2) + pow(U[i].Bz, 2) );
		W[i].P = W[i].P*(gamma-1);
	}
	
}
