#include <bits/stdc++.h> 
#include <iostream> 
#include <sys/stat.h> 
#include <sys/types.h> 

#include <fstream>
#include <stdlib.h>
#include <string>
#include "struct_vars.h"
#include "C_to_P.h"
#include "write_log.h"
using namespace std;

void log(struct U_VAR *U, struct B_VAR *B, int nx, int nghost, int test, double gamma){
    struct W_VAR W[nx];
    C_to_P(&U[0], &W[0], nx, gamma);

    // mode 0777 gives read/write/execute access to owner, group or other users
    mkdir("../log", 0777);	
	string dir, data_file;	
	switch (test){
	   case 1: // Brio Wu
         dir = "../log/Brio_Wu";
		 break;
	   case 2:
         dir = "../log/RJ2a";
		 break;
	   case 3:
         dir = "../log/Falle";
		 break;
	   case 4:
         dir = "../log/RJ4d";
		 break;
	   case 5:
         dir = "../log/Komissarov";
		 break;
	}	
	mkdir(&dir[0], 0777);	

	fstream file;
	
	data_file = dir +"/d.txt";
    file.open(&data_file[0], ios::out); //open file
    for(int i=nghost; i<nx-nghost; i++)
        file<<W[i].d<<"\n"; //write data to file
    file.close();	

	data_file = dir + "/ux.txt";
    file.open(&data_file[0], ios::out); //open file
    for(int i=nghost; i<nx-nghost; i++)
        file<<W[i].ux<<"\n"; //write data to file
    file.close();	
	
	data_file = dir + "/uy.txt";
    file.open(&data_file[0], ios::out); //open file
    for(int i=nghost; i<nx-nghost; i++)
        file<<W[i].uy<<"\n"; //write data to file
    file.close();	

	data_file = dir + "/uz.txt";
    file.open(&data_file[0], ios::out); //open file
    for(int i=nghost; i<nx-nghost; i++)
        file<<W[i].uz<<"\n"; //write data to file
    file.close();	
	
	data_file = dir + "/P.txt";
    file.open(&data_file[0], ios::out); //open file
    for(int i=nghost; i<nx-nghost; i++)
        file<<W[i].P<<"\n"; //write data to file
    file.close();	

	data_file = dir + "/Bx.txt";
    file.open(&data_file[0], ios::out); //open file
    for(int i=nghost; i<nx-nghost; i++)
        file<<W[i].Bx<<"\n"; //write data to file
    file.close();	

	data_file = dir + "/By.txt";
    file.open(&data_file[0], ios::out); //open file
    for(int i=nghost; i<nx-nghost; i++)
        file<<W[i].By<<"\n"; //write data to file
    file.close();	

	data_file = dir + "/Bz.txt";
    file.open(&data_file[0], ios::out); //open file
    for(int i=nghost; i<nx-nghost; i++)
        file<<W[i].Bz<<"\n"; //write data to file	
    file.close();	

}

/*if(!file){
    cerr << "Can't open file!\n";
    exit(1); // fail to open file
} */