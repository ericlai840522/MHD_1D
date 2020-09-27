#include<iostream>
#include "struct_vars.h"
#include "update_U.h"
using namespace std;

void update_U(struct U_VAR *U, struct U_VAR *Ut, struct U_VAR *F, int nx, double dx, double dt){
	
	int ip;
	for(int i=0;i<nx;i++){
		ip = i%nx + 1;
		Ut[i].d  = U[i].d  - dt/dx*( F[ip].d  - F[i].d );		
		Ut[i].Mx = U[i].Mx - dt/dx*( F[ip].Mx - F[i].Mx);		
		Ut[i].My = U[i].My - dt/dx*( F[ip].My - F[i].My);		
		Ut[i].Mz = U[i].Mz - dt/dx*( F[ip].Mz - F[i].Mz);		
		Ut[i].E  = U[i].E  - dt/dx*( F[ip].E  - F[i].E );		
		Ut[i].Bx = U[i].Bx - dt/dx*( F[ip].Bx - F[i].Bx);		
		Ut[i].By = U[i].By - dt/dx*( F[ip].By - F[i].By);		
		Ut[i].Bz = U[i].Bz - dt/dx*( F[ip].Bz - F[i].Bz);	
	}
}
