#include <cmath>
#include <iostream>
#include "Integrand_Classes.h"


using namespace std;

Plus_Distribution Curly_G_ia;

double My_Divergent(double* x){return x[0]*2+x[1]-x[2];}
double My_Finite(double* x){return x[0]*x[0];}


main(int argc, char* argv[]){        

    Curly_G_ia.set_D(&My_Divergent);
    Curly_G_ia.set_R(&My_Finite);
    
    double vars[4]={0.5,2,7,3};
    

    cout << "Curly_G_ia(0.5) = " << Curly_G_ia.eval_D(vars) << " , " << Curly_G_ia.eval_R(vars) << endl;


}                                      
