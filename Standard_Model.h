#ifndef _Standard_Model_H
#define _Standard_Model_H

#include <iostream>

using namespace std;

const int Nc = 3;

const double PISM = acos(-1.0);

const double Qu = 2.0/3.0;
const double Ql = -1.0;
const double Qw = 1.0;

///   Some Standard Model Constants   ////

  ////////////////////////////////////////////////////////
 ////                    Couplings                    ///
//////////////////////////////////////////////////////// 
 
 
const double cw = 80.376/91.1876;//80.385/91.1876;
const double sw = sqrt(1.0 - (cw*cw)) ;
const double e = sqrt(4.0*PISM/137.035999074);//0.30795376724436879;// Second setting to match MG5_aMC@NLO
const double g = sqrt (4.0*PISM*0.118); // At the mass of the Z


  //////////////////////////////////////////////////////
 ///                   b quark                      ///
//////////////////////////////////////////////////////

// const double mb =  4.6/;             //  In GeV/c^2
const double mb =  4.75;             //  In GeV/c^2
const double hc_b = -(1.0/2.0);     //  Hypercharge
const double ch_b = (1.0/3.0);      //  Charge


const double mz = 91.1876;//4.6/; //In GeV/c^2

const double cr = (1./3.)*(sw/cw);  /// Pr = (1/2)(1+Gamma5)
const double cl = ((-1./2)+(sw*sw/3))/(sw*cw); /// Pl = (1/2)(1-Gamm5)
const double me = 0.000510998928; /// Mass of the electron in GeV
const double mg = 0.0; /// Gluon Mass
const double mph = 0.0; /// Photon Mass
const double conversion = (0.000389379338)*1000000000000; // GeV^-2 to pBarn

void print_sm_parameters( void )
{
  
  cout << "Standard Model Parameters fixed at: \n" << "alpha = " << e*e/(4*PISM) << endl;
  
  cout << "cos(theta_w) = " << cw << endl;
  
  cout << "alpha_s = " << g*g/(4*PISM) << endl;
  
   
}



#endif