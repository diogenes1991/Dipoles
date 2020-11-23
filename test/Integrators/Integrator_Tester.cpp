#include "../../src/Integrators.h"

double my_f(double *k, size_t dim, void *params){return (k[0]*k[0]+k[1]);}
int my_g(const int *ndim, const double k[], const int *ncomp, double f[], void* params){f[0]=k[0]*k[0]+k[1];return 0;}

int main(){

    montecarlo_specs mc_sp;
    mc_sp.MaxEval = 100000;

    std::cout.precision(16);

    GSL_Integrator GIntegrator(my_f,2);
    GIntegrator.Integrate(&mc_sp,"Vegas");

    CUBA_Integrator CIntegrator(my_g,2);
    CIntegrator.Integrate(&mc_sp,"Vegas");


  return 0;
}