//////////////////////////////////////////////////////////////////////
//                                                                  //
//         T E S T   P R O G R A M                                  //
//                                                                  //
//////////////////////////////////////////////////////////////////////

#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <cmath>

// For C++ programs, the following needs to be included in order to use
// NLOX -  the only include needed to use the generated NLOX processes.
#include "nlox_olp.h"
#include "/home/diogenes1991/local/include/IR_Pole_NLOX.h"

// processconst.h is only included below such that the mass values from
// processconst.c can be accessed; only to convert the test PSP to BLHA 
// format. It isn't needed to actually use the generated NLOX processes
#include "processconst.h"

// Funtion to return all the color correlated matrix elements for a particular 
// coupling power combination



int main(int argc, char* argv[]) {
  using namespace std;
  ProcessConst pc;

  ////////////////////////////////////////////////////////////////////
  // PHASE SPACE POINT
  // (this is absolutely user/MC dependent, and only the form of the PSP
  // pp is relevant)

  // pp stores the phase-space point. It is an array of double precision
  // numbers of length 5 * (number of external particles).
  //
  // The elements of pp are
  //
  //    pp = [p1t, p1x, p1y, p1z, m1, p2t, p2x, ...]
  //
  // This is the format used by the BLHA2 standard. NLOX does not use
  // the mass component, so its value does not matter.
  //
  // There might be a few phase-space points listed here. Comment out 
  // all but the one you want to use.
  //
  // NOTE: The information in the mass entries and the last five 
  // entries of pp is currently not used by NLOX.

  // 2 to 2 Phase Space Point Generator Triggered with:
  // sqrt(s) = 500
  //      m1 = 0
  //      m2 = 0
  //      m3 = 0
  //      m4 = 0
  /// Standard NLOX format
  FourVector p1(250, 0 , 0 , 250);
  FourVector p2(250, 0 , 0 , -250);
  FourVector p3(250, -19.18429364851453 , 113.7028340258232 , 221.8189090490341);
  FourVector p4(250, 19.18429364851453 , -113.7028340258232 , -221.8189090490341);

  // Convert PSP to BLHA format.
  double pp[5*4] = { p1.p0, p1.p1, p1.p2, p1.p3, 0,
  p2.p0, p2.p1, p2.p2, p2.p3, 0,
  p3.p0, p3.p1, p3.p2, p3.p3, 0,
  p4.p0, p4.p1, p4.p2, p4.p3, 0 };
  int next = 4;

  ///
  ///
  ///
  /////////////////////////////////////////////////////////////////////
  
  
  ////////////////////////////////////////////////////////////////////
  ///    NLOX Initialization
  ///
  ///

  double acc2;
  double Alpha_s = 0.079577471545947667884441881686257181017229822870228;
  double Alpha_e = 0.079577471545947667884441881686257181017229822870228;
  double Pi = 3.1415926535897932384626433832795028841971693993751;
  
  int isub = 0;
  double mu = 500;
  int proc_id[4] = {2,-2,4,-4};
  
  NLOX_OLP_Start(NULL, NULL);
  
  bool verbose=1;
  
  ///
  ///
  ///
  ////////////////////////////////////////////////////////////////////////
  
  
  ///////////////////////////////////////////////////////////////
  ///  
  ///   We first compute the 3 possible color-correlated born matrix elements.
  ///
  
  char tt[]      = "tree_tree";
  char as2_ae0[] = "g2e0_g2e0";
  char as1_ae1[] = "g2e0_g0e2"; // We only compute one and multiply by 2 the result...
  char as0_ae2[] = "g0e2_g0e2";
  double B_as2_ae0[7],B_as1_ae1[7],B_as0_ae2[7];
  NLOX_OLP_EvalSubProcess_CC(&isub, tt, as2_ae0, pp, &next, &mu, B_as2_ae0, &acc2);
  NLOX_OLP_EvalSubProcess_CC(&isub, tt, as1_ae1, pp, &next, &mu, B_as1_ae1, &acc2);
  NLOX_OLP_EvalSubProcess_CC(&isub, tt, as0_ae2, pp, &next, &mu, B_as0_ae2, &acc2);
  
  for(int CCC=0;CCC<=6;CCC++)B_as1_ae1[CCC]*=2; // Here we multiply by 2 to account for the other combination
  
  if(verbose){
  cout << "alpha_s^2 alpha_e^0 LO" << endl;
  cout << B_as2_ae0[0] << endl;
  
  cout << "alpha_s^1 alpha_e^1 LO" << endl;
  cout << B_as1_ae1[0] << endl;
  
  cout << "alpha_s^0 alpha_e^2 LO" << endl;
  cout << B_as0_ae2[0] << endl;
  
  cout << endl;
  }
  
  
  ///
  ///
  ////////////////////////////////////////////////////////
  
  ////////////////////////////////////////////////////////
  /// 
  ///         Now we compute the NLO virtuals
  ///
  
  char tl[]        = "tree_loop";
  char as3_ae0[]   = "g2e0_g4e0";
  char as2_ae1_l[] = "g2e0_g2e2";
  char as2_ae1_r[] = "g0e2_g4e0";
  char as1_ae2_l[] = "g2e0_g0e4";
  char as1_ae2_r[] = "g0e2_g2e2";
  char as0_ae3[]   = "g0e2_g0e4";
  double V_as3_ae0[2],V_as2_ae1[2],V_as1_ae2[2],V_as0_ae3[2];
  double aux[2];
  
  
  NLOX_OLP_EvalSubProcess(&isub,tl,as3_ae0,pp,&next,&mu,V_as3_ae0,&acc2);
  
  if(verbose){
  cout << "alpha_s^3 alpha_e^0 Virtual NLO" <<endl;
  cout << as3_ae0 << endl;
  cout << "Double Pole = " << V_as3_ae0[0] << endl;
  cout << "Single Pole = " << V_as3_ae0[1] << endl;
  cout << endl;
  }
  
  
  NLOX_OLP_EvalSubProcess(&isub,tl,as2_ae1_l,pp,&next,&mu,V_as2_ae1,&acc2);
  
  if(verbose){
  cout << "alpha_s^2 alpha_e^1 Virtual NLO" <<endl;
  cout << as2_ae1_l << endl;
  cout << "Double Pole = " << V_as2_ae1[0] << endl;
  cout << "Single Pole = " << V_as2_ae1[1] << endl;
  cout << endl;
  }
  
  aux[0]=V_as2_ae1[0];aux[1]=V_as2_ae1[1];
  NLOX_OLP_EvalSubProcess(&isub,tl,as2_ae1_r,pp,&next,&mu,V_as2_ae1,&acc2);
  
  if(verbose){
  cout << as2_ae1_r << endl;
  cout << "Double Pole = " << V_as2_ae1[0] << endl;
  cout << "Single Pole = " << V_as2_ae1[1] << endl;
  cout << endl;
  }
  
  V_as2_ae1[0]+=aux[0];V_as2_ae1[1]+=aux[1];
  
  NLOX_OLP_EvalSubProcess(&isub,tl,as1_ae2_l,pp,&next,&mu,V_as1_ae2,&acc2);
  
  if(verbose){
  cout << "alpha_s^1 alpha_e^2 Virtual NLO" <<endl;
  cout << as1_ae2_l << endl;
  cout << "Double Pole = " << V_as1_ae2[0] << endl;
  cout << "Single Pole = " << V_as1_ae2[1] << endl;
  cout << endl;
  }
  
  aux[0]=V_as1_ae2[0];aux[1]=V_as1_ae2[1];
  NLOX_OLP_EvalSubProcess(&isub,tl,as1_ae2_r,pp,&next,&mu,V_as1_ae2,&acc2);
  
  if(verbose){
  cout << as1_ae2_r << endl;
  cout << "Double Pole = " << V_as1_ae2[0] << endl;
  cout << "Single Pole = " << V_as1_ae2[1] << endl;
  cout << endl;
  }
  
  V_as1_ae2[0]+=aux[0];V_as1_ae2[1]+=aux[1];
  
  NLOX_OLP_EvalSubProcess(&isub,tl,as0_ae3,pp,&next,&mu,V_as0_ae3,&acc2);
  
  if(verbose){
  cout << "alpha_s^0 alpha_e^3 Virtual NLO" <<endl;
  cout << as0_ae3 << endl;
  cout << "Double Pole = " << V_as0_ae3[0] << endl;
  cout << "Single Pole = " << V_as0_ae3[1] << endl;
  cout << endl;
  }
  
  ///
  ///
  ///  
  //////////////////////////////////////////////////////////////
  
  //////////////////////////////////////////////////////////////
  ///
  ///   Here we compute the 6 possible IR-Pole Combinations
  ///   Recall that Get_Process_EW_Poles_NEW  + (ae/2Pi)
  ///               Get_Process_QCD_Poles_NEW + (as/2Pi)  
  ///
  
  double R_as3_ae0[2],R_as2_ae1[2],R_as1_ae2[2],R_as0_ae3[2];
  
  
  Get_Process_QCD_Poles_NEW(next,proc_id,pp,mu,B_as2_ae0,R_as3_ae0);
  R_as3_ae0[0] *= (Alpha_s/(2*Pi));
  R_as3_ae0[1] *= (Alpha_s/(2*Pi));
  
  if(verbose){
  cout << "alpha_s^3 alpha_e^0 Real NLO" <<endl;
  cout << "QCD Dipoles on |QCD|^2 Born " << endl;
  cout << "Double Pole = " << R_as3_ae0[0] << endl;
  cout << "Single Pole = " << R_as3_ae0[1] << endl;
  cout << endl;
  }
  
  
  Get_Process_EW_Poles_NEW (next,proc_id,pp,mu,B_as2_ae0,R_as2_ae1);
  
  if(verbose){
  cout << "alpha_s^2 alpha_e^1 Real NLO" <<endl;
  cout << "EW Dipoles on |QCD|^2 Born " << endl;
  cout << "Double Pole = " << (Alpha_s/(2*Pi))*R_as2_ae1[0] << endl;
  cout << "Single Pole = " << (Alpha_s/(2*Pi))*R_as2_ae1[1] << endl;
  cout << endl;
  }
  
  aux[0]=(Alpha_e/(2*Pi))*R_as2_ae1[0];
  aux[1]=(Alpha_e/(2*Pi))*R_as2_ae1[1];
  Get_Process_QCD_Poles_NEW(next,proc_id,pp,mu,B_as1_ae1,R_as2_ae1);
  
  if(verbose){
  cout << "QCD Dipoles on 2Re(QCD*EW) Born " << endl;
  cout << "Double Pole = " << (Alpha_s/(2*Pi))*R_as2_ae1[0] << endl;
  cout << "Single Pole = " << (Alpha_s/(2*Pi))*R_as2_ae1[1] << endl;
  cout << endl;
  }
  
  R_as2_ae1[0]*=(Alpha_s/(2*Pi));
  R_as2_ae1[1]*=(Alpha_s/(2*Pi));
  
  R_as2_ae1[0]+=aux[0];
  R_as2_ae1[1]+=aux[1];
  
  Get_Process_EW_Poles_NEW (next,proc_id,pp,mu,B_as1_ae1,R_as1_ae2);
  
  if(verbose){
  cout << "alpha_s^1 alpha_e^2 Real NLO" <<endl;
  cout << "EW Dipoles on 2Re(QCD*EW) Born " << endl;
  cout << "Double Pole = " << (Alpha_s/(2*Pi))*R_as1_ae2[0] << endl;
  cout << "Single Pole = " << (Alpha_s/(2*Pi))*R_as1_ae2[1] << endl;
  cout << endl;
  }
  
  aux[0]=(Alpha_e/(2*Pi))*R_as1_ae2[0];
  aux[1]=(Alpha_e/(2*Pi))*R_as1_ae2[1];
  
  Get_Process_QCD_Poles_NEW(next,proc_id,pp,mu,B_as0_ae2,R_as1_ae2);
  
  if(verbose){
  cout << "QCD Dipoles on |EW|^2 Born " << endl;
  cout << "Double Pole = " << (Alpha_s/(2*Pi))*R_as1_ae2[0] << endl;
  cout << "Single Pole = " << (Alpha_s/(2*Pi))*R_as1_ae2[1] << endl;
  cout << endl;
  }
  
  R_as1_ae2[0]*=(Alpha_s/(2*Pi));
  R_as1_ae2[1]*=(Alpha_s/(2*Pi));
  
  R_as1_ae2[0]+=aux[0];
  R_as1_ae2[1]+=aux[1];
  
  Get_Process_EW_Poles_NEW (next,proc_id,pp,mu,B_as0_ae2,R_as0_ae3);
  R_as0_ae3[0]*=(Alpha_e/(2*Pi));
  R_as0_ae3[1]*=(Alpha_e/(2*Pi));
  
  if(verbose){
  cout << "alpha_s^0 alpha_e^3 Real NLO" <<endl;
  cout << "EW Dipoles on |EW|^2 Born " << endl;
  cout << "Double Pole = " << R_as0_ae3[0] << endl;
  cout << "Single Pole = " << R_as0_ae3[1] << endl;
  cout << endl;
  }
  
  ///
  ///
  ///
  /////////////////////////////////////////////////////////////////
  
  
  ////////////////////////////////////////////////////////////////
  ///
  ///                  Pole checking piece 
  ///
  
  cout << "alpha_s^3 alpha_e^0 NLO" <<endl;
  cout << "Virtual Double Pole = " << V_as3_ae0[0] << endl;
  cout << "Virtual Single Pole = " << V_as3_ae0[1] << endl;
  cout << "  Real  Double Pole = " << R_as3_ae0[0] << endl;
  cout << "  Real  Single Pole = " << R_as3_ae0[1] << endl;
  
  cout << "alpha_s^2 alpha_e^1 NLO" <<endl;
  cout << "Virtual Double Pole = " << V_as2_ae1[0] << endl;
  cout << "Virtual Single Pole = " << V_as2_ae1[1] << endl;
  cout << "  Real  Double Pole = " << R_as2_ae1[0] << endl;
  cout << "  Real  Single Pole = " << R_as2_ae1[1] << endl;
  
  cout << "alpha_s^1 alpha_e^2 NLO" <<endl;
  cout << "Virtual Double Pole = " << V_as1_ae2[0] << endl;
  cout << "Virtual Single Pole = " << V_as1_ae2[1] << endl;
  cout << "  Real  Double Pole = " << R_as1_ae2[0] << endl;
  cout << "  Real  Single Pole = " << R_as1_ae2[1] << endl;
  
  cout << "alpha_s^0 alpha_e^3 NLO" <<endl;
  cout << "Virtual Double Pole = " << V_as0_ae3[0] << endl;
  cout << "Virtual Single Pole = " << V_as0_ae3[1] << endl;
  cout << "  Real  Double Pole = " << R_as0_ae3[0] << endl;
  cout << "  Real  Single Pole = " << R_as0_ae3[1] << endl;
  
  ///
  ///
  ///
  //////////////////////////////////////////////////////////////////
  
  
  return 0;

}
