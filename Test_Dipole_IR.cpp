#include <iostream>
#include <fstream>
#include <cmath>
#include "Vector_Real.h"
#include "Phase_Spaces_Real.h"
#include "Dipole_Definitions.h"
#include "cuba.h"
#include <time.h>


#define NPAR 5

using namespace std;

double g_ab_no_dipole(double x, double s, double ma, double mb){
    
    double sbar = s - ma*ma - mb*mb;
    return (sbar*sbar)*(1.-x)/(16.*PI*PI*s);
    
}

int print_time(int time){
  int min = 60;
  int hour = 60*min;
  int day = 24*hour;
  
  float t_f;
  int t_i;
  
  t_f = float(time)/day;
  t_i = time/day;
  
  if ( t_i ) 
  {
    if ( t_i == 1 ) cout << t_i << " day ";
    else cout << t_i << " days ";
  }
    
  time = time - t_i*day;
  
  t_f = float(time)/hour;
  t_i = time/hour;
  
  if ( t_i )
  {
    if ( t_i == 1 ) cout << t_i << " hour ";
    else cout << t_i << " hours ";
  }
  
  time = time - t_i*hour;
  
  t_f = float(time)/min;
  t_i = time/min;
  
  if ( t_i )
    {
    if ( t_i == 1 ) cout << t_i << " minute ";
    else cout << t_i << " minutes ";
    }
  
  time = time - t_i*min;
  
  if ( time ) 
    {
    if ( time == 1 ) cout << time << " second";
    else cout << time << " seconds";
    }
  
  cout << endl;
  
  
}


/// Integration Targets to Check Total Volume of Phase Space ///

static int Integrand_PSV_Sub(const int *ndim, const double x[],const int *ncomp, double f[]){
    
    
    double Min[2] = {5,5};
    
    Vector Pfi[3];
    double Mfi[3] = {100,100,0};
    
    double fin[5] = {x[0],x[1],x[2],x[3],x[4]};
    
    double sqrts = 1000;
    
    Vector P(sqrts,0,0,0);
    
    double J1=1.;

    Recursive_PSP(P,3,Pfi,Mfi,fin,J1);
    
    double xab = 1.0 - ((Pfi[2]*P)/(P*P-Min[0]*Min[0]-Min[1]*Min[1]));
    
//     double xmin = 2*Min[0]*Min[1]/(P*P-Min[0]*Min[0]-Min[1]*Min[1]);
//     double xminalt = ((Mfi[0]+Mfi[1])*(Mfi[0]+Mfi[1]) - (Min[0]*Min[0]+Min[1]*Min[1]))/(P*P-Min[0]*Min[0]-Min[1]*Min[1]);
    
//     if (xmin<xminalt) xmin = xminalt;
    
    
    f[0] = J1;
    
    
}

static int Integrand_PSV_Plu(const int *ndim, const double x[],const int *ncomp, double f[]){
    
    double sqrts = 1000;
    
    double Min[2] = {5,5};
    
    double sbar = sqrts*sqrts - Min[0]*Min[0] - Min[1]*Min[1];
    double xab = (2*Min[0]*Min[1]/sbar)+x[0]*(1.-(2*Min[0]*Min[1]/sbar));

    double sqrtshat = sbar*xab + Min[0]*Min[0] + Min[1]*Min[1];
           sqrtshat = sqrt(sqrtshat);
    Vector Pabtil(sqrtshat,0,0,0);
    
    Vector Pfi[2];
    double Mfi[2]={100,100};
    
    if(sqrtshat>(Mfi[0]+Mfi[1])){
    
    double ficonf[2] = {x[1],x[2]};
    double Jfi = (1.-(2*Min[0]*Min[1]/sbar));
    
    Recursive_PSP(Pabtil,2,Pfi,Mfi,ficonf,Jfi);
    
    f[0] = Jfi*g_ab_no_dipole(xab,sqrts*sqrts,Min[0],Min[1]);
    }
    
    else {
        f[0] = 0.;
    }
        
    
    
    
}


/// Integration Targets to Check Analythical Integration of Dipole Functions                  ///
/// We employ a Dummy Matrix element set at x^2 to cancel any singularity upon regularization ///


static int Integrand_Dipole_Sub(const int *ndim, const double x[],const int *ncomp, double f[]){
    
    
    double Min[2] = {5,5};
    
    Vector Pfi[3];
    double Mfi[3] = {100,100,0};
    
    double fin[5] = {x[0],x[1],x[2],x[3],x[4]};
    
    double sqrts = 1000;
    
    Vector P(sqrts,0,0,0);
    
    double J1=1.;

    Recursive_PSP(P,3,Pfi,Mfi,fin,J1);
    
    double xab = 1.0 - ((Pfi[2]*P)/(P*P-Min[0]*Min[0]-Min[1]*Min[1]));
    
//     double xmin = 2*Min[0]*Min[1]/(P*P-Min[0]*Min[0]-Min[1]*Min[1]);
//     double xminalt = ((Mfi[0]+Mfi[1])*(Mfi[0]+Mfi[1]) - (Min[0]*Min[0]+Min[1]*Min[1]))/(P*P-Min[0]*Min[0]-Min[1]*Min[1]);
    
//     if (xmin<xminalt) xmin = xminalt;
    
    
    f[0] = J1;
    
    
}

static int Integrand_Dipole_Plu(const int *ndim, const double x[],const int *ncomp, double f[]){
    
    double sqrts = 1000;
    
    double Min[2] = {5,5};
    
    double sbar = sqrts*sqrts - Min[0]*Min[0] - Min[1]*Min[1];
    double xab = (2*Min[0]*Min[1]/sbar)+x[0]*(1.-(2*Min[0]*Min[1]/sbar));

    double sqrtshat = sbar*xab + Min[0]*Min[0] + Min[1]*Min[1];
           sqrtshat = sqrt(sqrtshat);
    Vector Pabtil(sqrtshat,0,0,0);
    
    Vector Pfi[2];
    double Mfi[2]={100,100};
    
    if(sqrtshat>(Mfi[0]+Mfi[1])){
    
    double ficonf[2] = {x[1],x[2]};
    double Jfi = (1.-(2*Min[0]*Min[1]/sbar));
    
    Recursive_PSP(Pabtil,2,Pfi,Mfi,ficonf,Jfi);
    
    f[0] = Jfi*g_ab_no_dipole(xab,sqrts*sqrts,Min[0],Min[1]);
    }
    
    else {
        f[0] = 0.;
    }
        
    
    
    
}

static int Integrand_Dipole_End(const int *ndim, const double x[],const int *ncomp, double f[]){
    
    double sqrts = 1000;
    
    double Min[2] = {5,5};
    
    double sbar = sqrts*sqrts - Min[0]*Min[0] - Min[1]*Min[1];
    double xab = (2*Min[0]*Min[1]/sbar)+x[0]*(1.-(2*Min[0]*Min[1]/sbar));

    double sqrtshat = sbar*xab + Min[0]*Min[0] + Min[1]*Min[1];
           sqrtshat = sqrt(sqrtshat);
    Vector Pabtil(sqrtshat,0,0,0);
    
    Vector Pfi[2];
    double Mfi[2]={100,100};
    
    if(sqrtshat>(Mfi[0]+Mfi[1])){
    
    double ficonf[2] = {x[1],x[2]};
    double Jfi = (1.-(2*Min[0]*Min[1]/sbar));
    
    Recursive_PSP(Pabtil,2,Pfi,Mfi,ficonf,Jfi);
    
    f[0] = Jfi*g_ab_no_dipole(xab,sqrts*sqrts,Min[0],Min[1]);
    }
    
    else {
        f[0] = 0.;
    }
        
    
    
    
}




/// Integartion Functions to perform the Integrals


void Integrate_Sub(){
    const int ndim      = 5;
    double userdata[1];
    const int ncomp     = 1;
    const int nvec      = 1;
    const double epsrel = 0;
    const double epsabs = 0;
    const int flags     = 1; // 1 = verbose
    const int seed      = 1;
    const int mineval   = 10000;
    const int maxeval   = 500000000;
    
    const int nstart    = 5000000;
    const int nincrease = 750000;
    const int nbatch    = 1000;
    const int gridno    = 0;
    const char* statefile = NULL;
    const int spin        = -1;
 
    int neval, fail;
    double integral[nvec], error[nvec], prob[nvec];
    
  
    int t;
    t = time(NULL);
//     cout << "Born Integration Initiated with seed = " << seed <<endl;
    
    Vegas(ndim, ncomp,
        Integrand_PSV_Sub, &userdata, nvec,
        epsrel, epsabs,
        flags, seed,
        mineval, maxeval,
        nstart, nincrease, nbatch,
        gridno, statefile, &spin,
        &neval, &fail,
        integral, error, prob);  

   cout << "Integral value over the n+1 Phase Space = (" << integral[0] << " +/- " <<  error[0]  << ") pb" << endl;

   print_time(time(NULL) - t);
   
   cout << "\n";
   
//   cout << "Ratio = " << integral[0]/0.000425041 << endl;
   
     
     
     
}

void Integrate_Plu(){
    const int ndim      = 3;
    double userdata[1];
    const int ncomp     = 1;
    const int nvec      = 1;
    const double epsrel = 0;
    const double epsabs = 0;
    const int flags     = 1; // 1 = verbose
    const int seed      = 1;
    const int mineval   = 10000;
    const int maxeval   = 500000000;
    
    const int nstart    = 5000000;
    const int nincrease = 750000;
    const int nbatch    = 1000;
    const int gridno    = 0;
    const char* statefile = NULL;
    const int spin        = -1;
 
    int neval, fail;
    double integral[nvec], error[nvec], prob[nvec];
    
  
    int t;
    t = time(NULL);
//     cout << "Born Integration Initiated with seed = " << seed <<endl;
    
    Vegas(ndim, ncomp,
        Integrand_PSV_Plu, &userdata, nvec,
        epsrel, epsabs,
        flags, seed,
        mineval, maxeval,
        nstart, nincrease, nbatch,
        gridno, statefile, &spin,
        &neval, &fail,
        integral, error, prob);  

   cout << "Integral value over the n(x) Phase Space = (" << integral[0] << " +/- " <<  error[0]  << ") pb" << endl;

   print_time(time(NULL) - t);
   
   cout << "\n";
   
//   cout << "Ratio = " << integral[0]/0.000425041 << endl;
   
     
     
     
}

void Integrate_End(){
    const int ndim      = 2;
    double userdata[1];
    const int ncomp     = 1;
    const int nvec      = 1;
    const double epsrel = 0;
    const double epsabs = 0;
    const int flags     = 1; // 1 = verbose
    const int seed      = 1;
    const int mineval   = 10000;
    const int maxeval   = 500000000;
    
    const int nstart    = 5000000;
    const int nincrease = 750000;
    const int nbatch    = 1000;
    const int gridno    = 0;
    const char* statefile = NULL;
    const int spin        = -1;
 
    int neval, fail;
    double integral[nvec], error[nvec], prob[nvec];
    
  
    int t;
    t = time(NULL);
//     cout << "Born Integration Initiated with seed = " << seed <<endl;
    
    Vegas(ndim, ncomp,
        Integrand_PSV_End, &userdata, nvec,
        epsrel, epsabs,
        flags, seed,
        mineval, maxeval,
        nstart, nincrease, nbatch,
        gridno, statefile, &spin,
        &neval, &fail,
        integral, error, prob);  

   cout << "Integral value over the n Phase Space = (" << integral[0] << " +/- " <<  error[0]  << ") pb" << endl;

   print_time(time(NULL) - t);
   
   cout << "\n";
   
//   cout << "Ratio = " << integral[0]/0.000425041 << endl;
   
     
     
     
}


main(int argc, char* argv[]){        



    Integrate_Sub();
    Integrate_Plu();
 

}     