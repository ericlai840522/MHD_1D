#ifndef INIT_PARAM_HEADER
#define INIT_PARAM_HEADER

#include <iostream>
#include <math.h>
#include "struct_vars.h"
using namespace std;
void init_param(int &test, struct U_VAR *U, struct B_VAR *B, int nx, double &tf, double &gamma){
    gamma = 5.0/3.0;
    switch (test){
        case 1: // Brio Wu
        // left hand side
        for(int i=0;i<nx/2;i++){ //the 1st to the (nx/2)th cells
            U[i].d = 1;
            U[i].Mx = 0;
            U[i].My = 0;
            U[i].Mz = 0;
            U[i].E = 2.28125;
            U[i].Bx = 0.75;
            U[i].By = 1;
            U[i].Bz = 0;
        }        
        // right hand side
        for(int i=nx/2; i<nx; i++){ //the (nx/2+1)th to the (nx)th cells
            U[i].d = 0.125;
            U[i].Mx = 0;
            U[i].My = 0;
            U[i].Mz = 0;
            U[i].E = 0.93125;
            U[i].Bx = 0.75;
            U[i].By = -1;
            U[i].Bz = 0;
        }
        tf = 0.1;
        break;
          
        case 2: // RJ2a
        // left hand side
        for(int i=0;i<nx/2;i++){ //the 1st to the (nx/2)th cells
            U[i].d = 1.08;
            U[i].Mx = 1.2*U[i].d;
            U[i].My = 0.01*U[i].d;
            U[i].Mz = 0.5*U[i].d;
            U[i].Bx = 2/sqrt(4*M_PI);
            U[i].By = 3.6/sqrt(4*M_PI);
            U[i].Bz = 2/sqrt(4*M_PI);
            U[i].E = 0.95/(gamma-1) + 0.5*( pow(U[i].Mx, 2) + pow(U[i].My, 2) + pow(U[i].Mz, 2) )/U[i].d + \
                0.5*( pow(U[i].Bx, 2) + pow(U[i].By, 2) + pow(U[i].Bz, 2) );
        }    

        // right hand side
        for(int i=nx/2; i<nx; i++){ //the (nx/2+1)th to the (nx)th cells
            U[i].d = 1;
            U[i].Mx = 0;
            U[i].My = 0;
            U[i].Mz = 0;
            U[i].Bx = 2/sqrt(4*M_PI);
            U[i].By = 4/sqrt(4*M_PI);
            U[i].Bz = 2/sqrt(4*M_PI);
            U[i].E = 1/(gamma-1) + 0.5*( pow(U[i].Mx, 2) + pow(U[i].My, 2) + pow(U[i].Mz, 2) )/U[i].d + \
                0.5*( pow(U[i].Bx, 2) + pow(U[i].By, 2) + pow(U[i].Bz, 2) );
         }    
        tf = 0.2;
        break;
        
        case 3: // Falle
        // left hand side
        for(int i=0;i<nx/2;i++){ //the 1st to the (nx/2)th cells
            U[i].d = 1.368;
            U[i].Mx = 0.269*U[i].d;
            U[i].My = 1*U[i].d;
            U[i].Mz = 0*U[i].d;
            U[i].Bx = 1;
            U[i].By = 0;
            U[i].Bz = 0;
            U[i].E = 3.88699;
        }    
        // right hand side
        for(int i=nx/2; i<nx; i++){ //the (nx/2+1)th to the (nx)th cells
            U[i].d = 1;
            U[i].Mx = 0;
            U[i].My = 0;
            U[i].Mz = 0;
            U[i].Bx = 1;
            U[i].By = 1;
            U[i].Bz = 0;
            U[i].E = 2.5;
        }    
        tf = 0.2;
        break;
         
        case 4: // RJ4d
        // left hand side
        for(int i=0;i<nx/2;i++){ //the 1st to the (nx/2)th cells
            U[i].d = 1;
            U[i].Mx = 0;
            U[i].My = 0;
            U[i].Mz = 0;
            U[i].Bx = 0.7;
            U[i].By = 0;
            U[i].Bz = 0;
            U[i].E = 1.745;
        }    
        // right hand side
        for(int i=nx/2; i<nx; i++){ //the (nx/2+1)th to the (nx)th cells
            U[i].d = 0.3;
            U[i].Mx = 0;
            U[i].My = 0;
            U[i].Mz = 1*U[i].d;
            U[i].Bx = 0.7;
            U[i].By = 1;
            U[i].Bz = 0;
            U[i].E = 1.195;
        }    
        tf = 0.15;
        break;
       
        case 5: // Komissarov
        // left hand side
        for(int i=0;i<nx/2;i++){ //the 1st to the (nx/2)th cells
            U[i].d = 1;
            U[i].Mx = 0;
            U[i].My = 0;
            U[i].Mz = 0;
            U[i].Bx = 1;
            U[i].By = 0;
            U[i].Bz = 0;
            U[i].E = 3.5;
        }    
        // right hand side
        for(int i=nx/2; i<nx; i++){ //the (nx/2+1)th to the (nx)th cells
            U[i].d = 0.2;
            U[i].Mx = 1.186/5;
            U[i].My = 2.967/5;
            U[i].Mz = 0;
            U[i].Bx = 1;
            U[i].By = 1.6405;
            U[i].Bz = 0;
            U[i].E = 3.0717;
        }    
        tf = 0.15; 
        break;
    }
    
    // face-centered magnetic field (left-faced)
    int im;
    for (int i=0; i<nx; i++){ 
        im = (i-1+nx)%nx; // im for i minus
        B[i].Bx = ( U[im].Bx + U[i].Bx )/2;
        B[i].By = ( U[im].By + U[i].By )/2;
        B[i].Bz = ( U[im].Bz + U[i].Bz )/2;
    }
    B[0].Bx = U[0].Bx;
    B[0].By = U[0].By;
    B[0].Bz = U[0].Bz;
}

#endif