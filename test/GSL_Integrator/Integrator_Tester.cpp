#include "../../src/Integrators.h"

struct MyArguments{
    size_t Dimension;
    size_t NCalls;
    size_t NIterations;
    void *Parameters;
    double (*Function)(double*,size_t,void*);

};

double my_f(double *k, size_t dim, void *params){return k[0]*k[0]+k[1];}

int main(){

    MyArguments IA;
    IA.Dimension = 2;
    IA.NCalls = 100000;
    IA.NIterations = 100;
    IA.Function = &my_f;

    GSL_Integrator<MyArguments> K(IA);
    K.Evaluate("Miser");



  return 0;
}