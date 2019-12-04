#include <iostream>
#include <fstream>
#include <cmath>
#include "Vector_Real.h"
#include "Phase_Spaces_Real.h"


main(int argc, char* argv[]){
    
const int NPar = atoi(argv[1]);
const double sqrts = atof(argv[2]);

Vector PList[NPar];
double MList[NPar];



for (int i=0;i<NPar;i++){MList[i]=atof(argv[i+3]);}


// Make sure P0 > Sum of all masses //
Vector P(sqrts,0,0,0);
double ran[3*NPar-4];

for (int i=0;i<(3*NPar-4);i++){ran[i]=aleatorio(0,1);}


double Jac=1;

Recursive_PSP(P,NPar,PList,MList,ran,Jac);

cout.precision(17);

Vector PTot;
for (int i=0;i<NPar;i++){
cout << "p" <<i+1 <<" = " << PList[i] << " m"<<i+1<<"^2 = " << PList[i]*PList[i] << endl;
PTot = PTot + PList[i];

}

cout << "P = " << PTot << endl;
cout << "Jacobian = " << Jac << endl;



}                                      
