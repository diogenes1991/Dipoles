#include <cmath>
#include <iostream>
#include "Integrand_Classes.h"


using namespace std;

Plus_Distribution Curly_G_ia;



main(int argc, char* argv[]){        


    cout << "Curly_G_ia(0.5) = " << Curly_G_ia.eval_D(0.5) << " , " << Curly_G_ia.eval_R(0.5) << endl;


}                                      
