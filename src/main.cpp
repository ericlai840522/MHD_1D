#include <iostream> 
#include <stdlib.h>
#include "struct_vars.h"
#include "init_param.h"
#include "cfl.h"
#include "LRstate.h"
#include "C_to_P.h"
#include "fluxes.h"
#include "HLLD.h"
#include "update_U.h"
#include "BC.h"
#include "write_log.h"
using namespace std;

int main(int argc, char *argv[]){ 
  int test = atoi(argv[1]); //convert char to int
  int nghost = 4; // number of ghost cell
  int nx = 800+2*nghost; // number of grid
  double Lx = 1; // domain length
  double dx = Lx/(nx-2*nghost); // grid size
  double t, tf, dt, gamma;
  int count = 0;
  
  struct U_VAR U[nx]; // cell-centered conserved variables
  struct U_VAR ULx[nx]; // left state of U
  struct U_VAR URx[nx]; // right state of U
  struct B_VAR B[nx]; // face-centered magnetic field
  struct W_VAR W[nx]; // Primitive variable
  struct U_VAR F[nx]; // face-centered flux
  
  init_param(test, &U[0], &B[0], nx, tf, gamma); //pass the memory location of U[0] to init
  dt = 0.8*cfl(&U[0], &W[0], nx, dx, gamma);
 
  struct U_VAR Uth[nx], Ut[nx];
  
  while(t<tf)
  {  
	  //step 1: compute the left and right state of U
	  LRstate(&ULx[0], &URx[0], &U[0], nx);
	  
	  //step 2: construct flux F
	  fluxes(&ULx[0], &URx[0], &F[0], &B[0], nx, gamma);	  
 	  
      //step 3: update U by a half time step
      update_U(&U[0],&Uth[0],&F[0],nx,dx,dt/2); // U^{t+dt/2}
	  
      //step 4: Apply boundary conditions	  
      BC(&Uth[0],&B[0],nghost,nx,test,gamma);


	  LRstate(&ULx[0], &URx[0], &Uth[0], nx);
	  fluxes(&ULx[0], &URx[0], &F[0], &B[0], nx, gamma); 	  
      update_U(&U[0],&Ut[0],&F[0],nx,dx,dt); // U^{t+dt}
	  for(int i=0; i<nx; i++)
		  U[i] = Ut[i];	
	  BC(&U[0],&B[0],nghost,nx,test,gamma);      

      t = t + dt;
      dt = 0.8*cfl(&U[0], &W[0], nx, dx, gamma);
      count = count + 1;	  
	  cout<<count<<" ";
  }
  
  log(&U[0], &B[0], nx, nghost, test, gamma);
  //for (int i=0;i<argc;i++)
  //	  cout<<argv[i]<<" ";

  return 0;
} 
