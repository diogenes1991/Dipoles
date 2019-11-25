#include <iostream>
#include <fstream>
#include <cmath>
#include "Vector_Real.h"
#include "Phase_Spaces_Real.h"

#define NPar 3

main(int argc, char* argv[]){        

Vector PList[NPar];
double MList[NPar];
// Make sure P0 > Sum of all masses //
Vector P(328,0,0,0);
double ran[3*NPar-4];

for (int i=0;i<NPar;i++){MList[i]=i+1;}
for (int i=0;i<(3*NPar-4);i++){ran[i]=aleatorio(0,1);}


Recursive_PSP(P,NPar,PList,MList,ran);

cout.precision(16);

Vector PTot;
for (int i=0;i<NPar;i++){
cout << "p" <<i+1 <<" = " << PList[i] << " m"<<i+1<<"^2 = " << PList[i]*PList[i] << endl;
PTot = PTot + PList[i];
}

cout << "P = " << PTot << endl;



}                                      
