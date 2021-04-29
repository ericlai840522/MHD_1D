#include <iostream>
#include "struct_vars.h"
#include "fluxes.h"
#include "HLLD.h"
using namespace std;

void fluxes(struct U_VAR *ULx, struct U_VAR *URx, struct U_VAR *F, struct B_VAR *B, int nx, double gamma){
    hlld(&ULx[0], &URx[0], &F[0], &B[0], nx, gamma);
}
