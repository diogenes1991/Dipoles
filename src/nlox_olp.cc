#include <cstring>
#include <cmath>
#include "nlox_olp.h"
#include "processconst.h"

Process* proc = NULL;

void NLOX_OLP_Start(char* fname, int* ierr) {
  proc = new Process();
}

void NLOX_OLP_EvalSubProcess_All(std::string SubProc, double* pp, int* next, double* mu, double* rval, double* acc) {
  int iv = *i;
  int nextv = *next;
  double muv = *mu;
  proc->evaluate_all(iv, pp, nextv, muv, rval, acc);
}

void NLOX_OLP_EvalSubProcess(std::string SubProc, char* type, char* cp, double* pp, int* next, double* mu, double* rval2, double* acc) {
  int iv = proc->AmpMap.at(SubProc);
  std::string types = type;
  std::string cps = cp;
  int nextv = *next;
  double muv = *mu;
  // Make sure only the Born is requested
  std::string ccLStrings[7] = {"ccL0","ccL1","ccL2","ccL3","ccL4","ccL5","ccL6"};
  proc->pc.set_param(ccLStrings[0],1.);
  for(int j=1; j<=nextv; j++) proc->pc.set_param(ccLStrings[j],0.);
  // Evaluate
  if(cps.find("as")==0) proc->evaluate_alpha(iv, types, cps, pp, nextv, muv, rval2, acc);
  else proc->evaluate(iv, types, cps, pp, nextv, muv, rval2, acc);
}

void NLOX_OLP_EvalSubProcess_CC(std::string SubProc, char* type, char* cp, double* pp, int* next, double* mu, double* rvalcc, double* acc) {
  int iv = *i;
  std::string types = type;
  std::string cps = cp;
  int nextv = *next;
  double muv = *mu;

  if ( types == "tree_loop" ) {
    std::cout << "Error:" << std::endl;
    std::cout << "Color correlated matrix elements are not generated for the specified intereference type \"" << types << "\"." << std::endl;
    abort();
  }

  // For n=next external particles, the return array rvalcc will be of the form 
  // {Born, <T1T2>, <T1T3>,...,<T1Tn>,<T2T3>,<T2T4>,...,<T2Tn>,...,<T(n-1)Tn>}.

  std::string ccLStrings[7] = {"ccL0","ccL1","ccL2","ccL3","ccL4","ccL5","ccL6"};
  double rval2[3]; 

  // Set and evaluate the Born
  proc->pc.set_param(ccLStrings[0],1.);
  for(int j=1; j<=nextv; j++) proc->pc.set_param(ccLStrings[j],0.);
  if(cps.find("as")==0) proc->evaluate_alpha(iv, types, cps, pp, nextv, muv, rval2, acc);
  else  proc->evaluate(iv, types, cps, pp, nextv, muv, rval2, acc);
  rvalcc[0] = rval2[2];
  
  // Set and evaluate <TiTj>
  int cccounter = 1;
  proc->pc.set_param(ccLStrings[0],0.);
  for(int j=1; j<=nextv; j++) {
    for(int k=j+1; k<=nextv; k++) {
      for(int l=0; l<=6; l++) {
        // Request <TjTk>=<TkTj>, k>j
        if(l==j||l==k) proc->pc.set_param(ccLStrings[l],1.);
        else proc->pc.set_param(ccLStrings[l],0.);
      }
      // Evaluate
      if(cps.find("as")==0) proc->evaluate_alpha(iv, types, cps, pp, nextv, muv, rval2, acc);
      else proc->evaluate(iv, types, cps, pp, nextv, muv, rval2, acc);
      rvalcc[cccounter] = rval2[2];
      cccounter += 1;
    }
  } 
}

void NLOX_OLP_SetMass(char* name, double* mass, double* width, int* ierr ){
    if(proc!=NULL){
        std::cout << "Warning: you are attempting to change a mass parameter after NLOX has been initialized,\n this change will not have an effect. Please set all mass parameters prior to the NLOX initialization call, or use NLOX_OLP_SetParameter to change at runtime.\n";
        *ierr = 0;
    }
    else{
    std::string mass_name(name);
    std::string width_name = mass_name;
    width_name.replace(0,1,"w");
    ProcessConst::set_param_in_file(mass_name,*mass,"nlox_parameters.par");
    ProcessConst::set_param_in_file(width_name,*width,"nlox_parameters.par");
    *ierr = 1;
    }   
}

void NLOX_OLP_SetParameter(char* para, double* re, double* im, int* ierr) {
  // TODO: check exact implementation standards of OLP. In principle para
  // is set by contract file, but we will default to using our own names.

  if (*im != 0.) { // No settable parameter is supported as complex in NLOX
    std::cout << "Warning: attempt to set unsupported imaginary part of parameter."
              << std::endl;
  }

  *ierr = proc->pc.set_param(std::string(para), *re);
  proc->update_mass(proc->pc.mass_name(para));
}
