#ifndef STRUCT_VARS_HEADER
#define STRUCT_VARS_HEADER

struct U_VAR{ //cell-centered
    double d; //density
    double Mx; //momentum x   
    double My; //momentum y
    double Mz; //momentum z
    double E; //total energy
    double Bx; //magnetic field x
    double By; //magnetic field y
    double Bz; //magnetic field z
}; 

struct B_VAR{ //face-centered
    double Bx; //magnetic field x
    double By; //magnetic field y
    double Bz; //magnetic field z
};

struct W_VAR{ //cell-centered
    double d; //density
    double ux; //velocity x   
    double uy; //velocity y
    double uz; //velocity z
    double P; //gas pressure
    double Bx; //magnetic field x
    double By; //magnetic field y
    double Bz; //magnetic field z
}; 

#endif