#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <string>
#include "struct_vars.h"
#include "C_to_P.h"
#include "fluxes.h"
#include "HLLD.h"
using namespace std;

void hlld(struct U_VAR *Ul, struct U_VAR *Ur, struct U_VAR *F, struct B_VAR *B, int nx, double gamma){
    
    string type;
    type = "Rusanov";
    type = "HLL";
    type = "HLLD";
    
    double bx[nx];
    for(int i=0; i<nx; i++){
        bx[i] = B[i].Bx;  
    }
    
    double spd[5][nx];
    
    struct U_VAR Fl[nx]; // left state of F
    struct U_VAR Fr[nx]; // right state of F
  
    struct U_VAR Ulst[nx]; // st for star in the paper Miyoshi & Kusano
    struct U_VAR Urst[nx]; // st for star in the paper Miyoshi & Kusano
    struct U_VAR Uldst[nx]; // dst for double star in the paper Miyoshi & Kusano
    struct U_VAR Urdst[nx]; // dst for double star in the paper Miyoshi & Kusano
  
    struct W_VAR Wl[nx];
    struct W_VAR Wr[nx];
    struct W_VAR Wlst[nx]; // st for star in the paper Miyoshi & Kusano
    struct W_VAR Wrst[nx]; // st for star in the paper Miyoshi & Kusano
    
    for(int i=0; i<nx; i++){
        Fl[i].d = 0;
        Fl[i].Mx = 0;
        Fl[i].My = 0;
        Fl[i].Mz = 0;
        Fl[i].E = 0;
        Fl[i].Bx = 0;
        Fl[i].By = 0;
        Fl[i].Bz = 0;
    }
    
    for(int i=0; i<nx; i++){
        Fr[i]    = Fl[i];
        Ulst[i]  = Fl[i];
        Urst[i]  = Fl[i];
        Uldst[i] = Fl[i];
        Urdst[i] = Fl[i];
    }
     
    C_to_P(&Ul[0], &Wl[0], nx, gamma);
    C_to_P(&Ur[0], &Wr[0], nx, gamma);
  
    // Step 1:  Compute left & right wave speeds 
    double pbl[nx], pbr[nx], gpl[nx], gpr[nx], gpbl[nx], gpbr[nx];
    for(int i=0; i<nx; i++){
        pbl[i] = 0.5*(pow(bx[i], 2) + pow(Wl[i].By, 2) + pow(Wl[i].Bz, 2));
        pbr[i] = 0.5*(pow(bx[i], 2) + pow(Wr[i].By, 2) + pow(Wr[i].Bz, 2));      
        gpl[i]  = gamma*Wl[i].P;
        gpr[i]  = gamma*Wr[i].P;
        gpbl[i] = gpl[i] + 2.0*pbl[i];
        gpbr[i] = gpr[i] + 2.0*pbr[i];      
    }
     
    double cfl[nx], cfr[nx], cfmax[nx], smax[nx];
    for(int i=0; i<nx; i++){
        cfl[i] = sqrt( (gpbl[i] + sqrt( pow(gpbl[i], 2) - 4.0*gpl[i]*pow(bx[i], 2)))/(2.0*Wl[i].d) ); // Miyoshi & Kusano, eqn. (3)
        cfr[i] = sqrt( (gpbr[i] + sqrt( pow(gpbr[i], 2) - 4.0*gpr[i]*pow(bx[i], 2)))/(2.0*Wr[i].d) ); 
        cfmax[i] = (cfl[i] > cfr[i]) ? cfl[i] : cfr[i];
        smax[i]  =  max( fabs(Wl[i].ux) + cfl[i], fabs(Wr[i].ux) + cfr[i] );      
    }
   
    // This approach of determining spd can lead to noises in certain area
    for(int i=0; i<nx; i++){
        spd[0][i] = min(Wl[i].ux, Wr[i].ux) - max(cfl[i], cfr[i]);
        spd[4][i] = max(Wl[i].ux, Wr[i].ux) + max(cfl[i], cfr[i]);
    }
  
    // eqn.(10.56) in Toro, 
    // Riemann Solvers and Numerical Methods for Fluid Dynamics
    for(int i=0; i<nx; i++){      
        spd[4][i] = max( fabs(Wl[i].ux) + cfl[i], fabs(Wr[i].ux) + cfr[i] );
        spd[0][i] = -spd[4][i];
    }
    
    // Step 2: Compute L/R fluxes
    double ptl[nx], ptr[nx];
    for(int i=0; i<nx; i++){
        ptl[i] = Wl[i].P + pbl[i]; // total pressure
        ptr[i] = Wr[i].P + pbr[i];
    }  
  
    for(int i=0; i<nx; i++){
        Fl[i].d  = Ul[i].Mx;
        Fl[i].Mx = Ul[i].Mx*Wl[i].ux + ptl[i] - pow(bx[i], 2);
        Fl[i].My = Ul[i].Mx*Wl[i].uy - bx[i]*Ul[i].By;
        Fl[i].Mz = Ul[i].Mx*Wl[i].uz - bx[i]*Ul[i].Bz;
        Fl[i].E  = (Ul[i].E + ptl[i])*Wl[i].ux - bx[i]*(Wl[i].ux*bx[i] + Wl[i].uy*Ul[i].By + Wl[i].uz*Ul[i].Bz);
        Fl[i].Bx = 0;
        Fl[i].By = Wl[i].ux*Ul[i].By - Wl[i].uy*bx[i];
        Fl[i].Bz = Wl[i].ux*Ul[i].Bz - Wl[i].uz*bx[i];
  
        Fr[i].d  = Ur[i].Mx;
        Fr[i].Mx = Ur[i].Mx*Wr[i].ux + ptr[i] - pow(bx[i], 2);
        Fr[i].My = Ur[i].Mx*Wr[i].uy - bx[i]*Ur[i].By;
        Fr[i].Mz = Ur[i].Mx*Wr[i].uz - bx[i]*Ur[i].Bz;
        Fr[i].E  = (Ur[i].E + ptr[i])*Wr[i].ux - bx[i]*(Wr[i].ux*bx[i] + Wr[i].uy*Ur[i].By + Wr[i].uz*Ur[i].Bz);
        Fr[i].Bx = 0;
        Fr[i].By = Wr[i].ux*Ur[i].By - Wr[i].uy*bx[i];
        Fr[i].Bz = Wr[i].ux*Ur[i].Bz - Wr[i].uz*bx[i];
    }
  
    // Step 3: Compute middle and Alfven wave speeds  
    // eqn (38) of Miyoshi & Kusano
    double sdl[nx], sdr[nx], sdml[nx], sdmr[nx], ptst[nx];
    for(int i=0; i<nx; i++){  
        sdl[i] = spd[0][i] - Wl[i].ux;
        sdr[i] = spd[4][i] - Wr[i].ux;
        spd[2][i] = (sdr[i]*Wr[i].d*Wr[i].ux - sdl[i]*Wl[i].d*Wl[i].ux - ptr[i] + ptl[i])/(sdr[i]*Wr[i].d - sdl[i]*Wl[i].d); // Sm
        sdml[i] = spd[0][i] - spd[2][i];
        sdmr[i] = spd[4][i] - spd[2][i];
    }
  
    // Step 4: Compute intermediate states Ul* 
    // eqn (43) of Miyoshi & Kusano
    for(int i=0; i<nx; i++){  
        Ulst[i].d = Ul[i].d*sdl[i]/sdml[i];
        Urst[i].d = Ur[i].d*sdr[i]/sdmr[i];
        // eqn (51) of Miyoshi & Kusano
        spd[1][i] = spd[2][i] - fabs(bx[i])/sqrt(Ulst[i].d); //fabs is defined in "stdlib.h"
        spd[3][i] = spd[2][i] + fabs(bx[i])/sqrt(Urst[i].d);
        ptst[i] = ptl[i] + Ul[i].d*sdl[i]*(sdl[i]-sdml[i]); // eqn(23)
    }
    
    // eqn (39) of Miyoshi & Kusano
    for(int i=0; i<nx; i++)
        Ulst[i].Mx = Ulst[i].d*spd[2][i];
  
    double small_num = 1.0e-100;
    double tmp;
    for(int i=0; i<nx; i++){  
        if ( fabs( Ul[i].d*sdl[i]*sdml[i]-pow(bx[i],2) ) < small_num*ptst[i] ){
            Ulst[i].My = Ulst[i].d*Wl[i].uy;
            Ulst[i].Mz = Ulst[i].d*Wl[i].uz;
            Ulst[i].By = Ul[i].By;
            Ulst[i].Bz = Ul[i].Bz;
        }
        else{ 
            // eqns (44) and (46) of M&K
            tmp = bx[i]*(sdl[i]-sdml[i])/(Ul[i].d*sdl[i]*sdml[i] - pow(bx[i],2) + small_num);
            Ulst[i].My = Ulst[i].d*( Wl[i].uy - Ul[i].By*tmp );
            Ulst[i].Mz = Ulst[i].d*( Wl[i].uz - Ul[i].Bz*tmp );
            // eqns (45) and (47) of M&K
            tmp = ( Ul[i].d*pow(sdl[i],2) - pow(bx[i],2) )/(Ul[i].d*sdl[i]*sdml[i] - pow(bx[i],2) + small_num);
            Ulst[i].By = Ul[i].By*tmp;
            Ulst[i].Bz = Ul[i].Bz*tmp;
        }
    }
    
    // eqns (48) of M&K
    double vbstl[nx];
    for(int i=0; i<nx; i++){
        vbstl[i] = (Ulst[i].Mx*bx[i] + Ulst[i].My*Ulst[i].By + Ulst[i].Mz*Ulst[i].Bz)/Ulst[i].d;
        Ulst[i].E = ( sdl[i]*Ul[i].E - ptl[i]*Wl[i].ux + ptst[i]*spd[2][i] + bx[i]*(Wl[i].ux*bx[i] + Wl[i].uy*Ul[i].By + Wl[i].uz*Ul[i].Bz - vbstl[i]) )/sdml[i];
    }
    C_to_P(&Ulst[0], &Wlst[0], nx, gamma);
  
    // Ur*
    // eqn (39) of Miyoshi & Kusano
    for(int i=0; i<nx; i++)
        Urst[i].Mx = Urst[i].d*spd[2][i];
      
    for(int i=0; i<nx; i++){  
        if ( fabs( Ur[i].d*sdr[i]*sdmr[i]-pow(bx[i],2) ) < small_num*ptst[i] ){
            Urst[i].My = Urst[i].d*Wr[i].uy;
            Urst[i].Mz = Urst[i].d*Wr[i].uz;
            Urst[i].By = Ur[i].By;
            Urst[i].Bz = Ur[i].Bz;
        }
        else{ 
            // eqns (44) and (46) of M&K
            tmp = bx[i]*(sdr[i]-sdmr[i])/(Ur[i].d*sdr[i]*sdmr[i] - pow(bx[i],2) + small_num);
            Urst[i].My = Urst[i].d*( Wr[i].uy - Ur[i].By*tmp );
            Urst[i].Mz = Urst[i].d*( Wr[i].uz - Ur[i].Bz*tmp );
            // eqns (45) and (47) of M&K
            tmp = ( Ur[i].d*pow(sdr[i],2) - pow(bx[i],2) )/(Ur[i].d*sdr[i]*sdmr[i] - pow(bx[i],2) + small_num);
            Urst[i].By = Ur[i].By*tmp;
            Urst[i].Bz = Ur[i].Bz*tmp;
        }
    }
    
    // eqns (48) of M&K
    double vbstr[nx];
    for(int i=0; i<nx; i++){
        vbstr[i] = (Urst[i].Mx*bx[i] + Urst[i].My*Urst[i].By + Urst[i].Mz*Urst[i].Bz)/Urst[i].d;
        Urst[i].E = ( sdr[i]*Ur[i].E - ptr[i]*Wr[i].ux + ptst[i]*spd[2][i] + bx[i]*(Wr[i].ux*bx[i] + Wr[i].uy*Ur[i].By + Wr[i].uz*Ur[i].Bz - vbstr[i]) )/sdmr[i];
    }
    C_to_P(&Urst[0], &Wrst[0], nx, gamma);
  
    // Ul** and Ur** - if Bx is zero, same as *-states 
    double invsumd;
    int Bxsign;
    for(int i=0; i<nx; i++){
        if(0.5*pow(bx[i],2) < small_num*ptst[i]){
            Uldst[i] = Ulst[i]; // all of the eight variables in the struct are copied
            Urdst[i] = Urst[i]; // all of the eight variables in the struct are copied
        }
        else{
            invsumd = 1.0/(sqrt(Ulst[i].d) + sqrt(Urst[i].d));
     
            if(bx[i] > 0.0) 
                Bxsign =  1.0;
            else if(bx[i] == 0.0)
                Bxsign = 0.0;
            else
                Bxsign = -1.0;
          
            Uldst[i].d = Ulst[i].d;
            Urdst[i].d = Urst[i].d;
            Uldst[i].Mx = Ulst[i].Mx;
            Urdst[i].Mx = Urst[i].Mx;
            // eqn (59) of M&K
            tmp = invsumd*(sqrt(Ulst[i].d)*Wlst[i].uy + sqrt(Urst[i].d)*Wrst[i].uy + Bxsign*(Urst[i].By-Ulst[i].By));
            Uldst[i].My = Uldst[i].d*tmp;
            Urdst[i].My = Urdst[i].d*tmp;
            // eqn (60) of M&K
            tmp = invsumd*(sqrt(Ulst[i].d)*Wlst[i].uz + sqrt(Urst[i].d)*Wrst[i].uz + Bxsign*(Urst[i].Bz-Ulst[i].Bz));
            Uldst[i].Mz = Uldst[i].d*tmp;
            Urdst[i].Mz = Urdst[i].d*tmp;
            // eqn (61) of M&K
            tmp = invsumd*(sqrt(Ulst[i].d)*Urst[i].By + sqrt(Urst[i].d)*Ulst[i].By + \
                            Bxsign*sqrt(Ulst[i].d)*sqrt(Urst[i].d)*(Wrst[i].uy-Wlst[i].uy));
            Uldst[i].By = tmp;
            Urdst[i].By = tmp;
            // eqn (62) of M&K
            tmp = invsumd*(sqrt(Ulst[i].d)*Urst[i].Bz + sqrt(Urst[i].d)*Ulst[i].Bz + \
                            Bxsign*sqrt(Ulst[i].d)*sqrt(Urst[i].d)*(Wrst[i].uz-Wlst[i].uz));
            Uldst[i].Bz = tmp;
            Urdst[i].Bz = tmp;
            // eqn (63) of M&K
            tmp = spd[2][i]*bx[i] + (Uldst[i].My*Uldst[i].By + Uldst[i].Mz*Uldst[i].Bz)/Uldst[i].d;
            Uldst[i].E = Ulst[i].E - sqrt(Ulst[i].d)*Bxsign*(vbstl[i] - tmp);
            Urdst[i].E = Urst[i].E + sqrt(Urst[i].d)*Bxsign*(vbstr[i] - tmp);
       } 
    }
    
    if(type == "Rusanov"){
        // Rusanov solver 
        for(int i=0; i<nx; i++){
            F[i].d  = 0.5*( Fl[i].d  + Fr[i].d  - smax[i]*(Ur[i].d  - Ul[i].d ) );
            F[i].Mx = 0.5*( Fl[i].Mx + Fr[i].Mx - smax[i]*(Ur[i].Mx - Ul[i].Mx) );
            F[i].My = 0.5*( Fl[i].My + Fr[i].My - smax[i]*(Ur[i].My - Ul[i].My) );
            F[i].Mz = 0.5*( Fl[i].Mz + Fr[i].Mz - smax[i]*(Ur[i].Mz - Ul[i].Mz) );
            F[i].E  = 0.5*( Fl[i].E  + Fr[i].E  - smax[i]*(Ur[i].E  - Ul[i].E ) );
            F[i].By = 0.5*( Fl[i].By + Fr[i].By - smax[i]*(Ur[i].By - Ul[i].By) );
            F[i].Bz = 0.5*( Fl[i].Bz + Fr[i].Bz - smax[i]*(Ur[i].Bz - Ul[i].Bz) );
        }
    }
  
    else if(type == "HLL"){ 
        // HLL solver 
        // Return upwind flux if flow is supersonic
        for(int i=0; i<nx; i++){
            if(spd[0][i] > 0.0)
                F[i] = Fl[i];
            else if(spd[4][i] < 0.0)
                F[i] = Fr[i];
            else if (spd[0][i] <= 0.0 && spd[4][i] >= 0.0){
                F[i].d  = (spd[4][i]*Fl[i].d  - spd[0][i]*Fr[i].d  + spd[0][i]*spd[4][i]*(Ur[i].d  - Ul[i].d ))/(spd[4][i] - spd[0][i]);
                F[i].Mx = (spd[4][i]*Fl[i].Mx - spd[0][i]*Fr[i].Mx + spd[0][i]*spd[4][i]*(Ur[i].Mx - Ul[i].Mx))/(spd[4][i] - spd[0][i]);
                F[i].My = (spd[4][i]*Fl[i].My - spd[0][i]*Fr[i].My + spd[0][i]*spd[4][i]*(Ur[i].My - Ul[i].My))/(spd[4][i] - spd[0][i]);
                F[i].Mz = (spd[4][i]*Fl[i].Mz - spd[0][i]*Fr[i].Mz + spd[0][i]*spd[4][i]*(Ur[i].Mz - Ul[i].Mz))/(spd[4][i] - spd[0][i]);
                F[i].E  = (spd[4][i]*Fl[i].E  - spd[0][i]*Fr[i].E  + spd[0][i]*spd[4][i]*(Ur[i].E  - Ul[i].E ))/(spd[4][i] - spd[0][i]);
                F[i].By = (spd[4][i]*Fl[i].By - spd[0][i]*Fr[i].By + spd[0][i]*spd[4][i]*(Ur[i].By - Ul[i].By))/(spd[4][i] - spd[0][i]);
                F[i].Bz = (spd[4][i]*Fl[i].Bz - spd[0][i]*Fr[i].Bz + spd[0][i]*spd[4][i]*(Ur[i].Bz - Ul[i].Bz))/(spd[4][i] - spd[0][i]);
            }
        }
    }
  
    else if(type == "HLLD"){  
        // HLLD solver 
        for(int i=0; i<nx; i++){
            if(spd[0][i] > 0.0){
                F[i].d  = Fl[i].d ;
                F[i].Mx = Fl[i].Mx;
                F[i].My = Fl[i].My;
                F[i].Mz = Fl[i].Mz;
                F[i].E  = Fl[i].E ;
                F[i].By = Fl[i].By;
                F[i].Bz = Fl[i].Bz;
            }         
            if(spd[4][i] < 0.0){
                F[i].d  = Fr[i].d ;
                F[i].Mx = Fr[i].Mx;
                F[i].My = Fr[i].My;
                F[i].Mz = Fr[i].Mz;
                F[i].E  = Fr[i].E ;
                F[i].By = Fr[i].By;
                F[i].Bz = Fr[i].Bz;         
            }
           
            if(spd[0][i] <= 0.0 && spd[1][i] >= 0.0){ // if SL = 0.0 ,then case 2 is reduced to case 1
                F[i].d  = Fl[i].d   + spd[0][i]*(Ulst[i].d  - Ul[i].d);
                F[i].Mx = Fl[i].Mx  + spd[0][i]*(Ulst[i].Mx - Ul[i].Mx);
                F[i].My = Fl[i].My  + spd[0][i]*(Ulst[i].My - Ul[i].My);
                F[i].Mz = Fl[i].Mz  + spd[0][i]*(Ulst[i].Mz - Ul[i].Mz);
                F[i].E  = Fl[i].E   + spd[0][i]*(Ulst[i].E  - Ul[i].E);
                F[i].By = Fl[i].By  + spd[0][i]*(Ulst[i].By - Ul[i].By);
                F[i].Bz = Fl[i].Bz  + spd[0][i]*(Ulst[i].Bz - Ul[i].Bz);
            }
            else if(spd[2][i] >= 0.0){ // case 3
                tmp = spd[1][i] - spd[0][i];
                F[i].d  = Fl[i].d  - spd[0][i]*Ul[i].d  - tmp*Ulst[i].d  + spd[1][i]*Uldst[i].d;
                F[i].Mx = Fl[i].Mx - spd[0][i]*Ul[i].Mx - tmp*Ulst[i].Mx + spd[1][i]*Uldst[i].Mx;
                F[i].My = Fl[i].My - spd[0][i]*Ul[i].My - tmp*Ulst[i].My + spd[1][i]*Uldst[i].My;
                F[i].Mz = Fl[i].Mz - spd[0][i]*Ul[i].Mz - tmp*Ulst[i].Mz + spd[1][i]*Uldst[i].Mz;
                F[i].E  = Fl[i].E  - spd[0][i]*Ul[i].E  - tmp*Ulst[i].E  + spd[1][i]*Uldst[i].E;
                F[i].By = Fl[i].By - spd[0][i]*Ul[i].By - tmp*Ulst[i].By + spd[1][i]*Uldst[i].By;
                F[i].Bz = Fl[i].Bz - spd[0][i]*Ul[i].Bz - tmp*Ulst[i].Bz + spd[1][i]*Uldst[i].Bz;
            }
            else if(spd[3][i] >= 0.0){ // case 4
                tmp = spd[3][i] - spd[4][i];
                F[i].d  = Fr[i].d  - spd[4][i]*Ur[i].d  - tmp*Urst[i].d  + spd[3][i]*Urdst[i].d;
                F[i].Mx = Fr[i].Mx - spd[4][i]*Ur[i].Mx - tmp*Urst[i].Mx + spd[3][i]*Urdst[i].Mx;
                F[i].My = Fr[i].My - spd[4][i]*Ur[i].My - tmp*Urst[i].My + spd[3][i]*Urdst[i].My;
                F[i].Mz = Fr[i].Mz - spd[4][i]*Ur[i].Mz - tmp*Urst[i].Mz + spd[3][i]*Urdst[i].Mz;
                F[i].E  = Fr[i].E  - spd[4][i]*Ur[i].E  - tmp*Urst[i].E  + spd[3][i]*Urdst[i].E;
                F[i].By = Fr[i].By - spd[4][i]*Ur[i].By - tmp*Urst[i].By + spd[3][i]*Urdst[i].By;
                F[i].Bz = Fr[i].Bz - spd[4][i]*Ur[i].Bz - tmp*Urst[i].Bz + spd[3][i]*Urdst[i].Bz;
            }
            else if(spd[4][i] >= 0.0){ // case 5
                F[i].d  = Fr[i].d  + spd[4][i]*( Urst[i].d  - Ur[i].d);
                F[i].Mx = Fr[i].Mx + spd[4][i]*( Urst[i].Mx - Ur[i].Mx);
                F[i].My = Fr[i].My + spd[4][i]*( Urst[i].My - Ur[i].My);
                F[i].Mz = Fr[i].Mz + spd[4][i]*( Urst[i].Mz - Ur[i].Mz);
                F[i].E  = Fr[i].E  + spd[4][i]*( Urst[i].E  - Ur[i].E);
                F[i].By = Fr[i].By + spd[4][i]*( Urst[i].By - Ur[i].By);
                F[i].Bz = Fr[i].Bz + spd[4][i]*( Urst[i].Bz - Ur[i].Bz);
            }
        }
    }
  
}