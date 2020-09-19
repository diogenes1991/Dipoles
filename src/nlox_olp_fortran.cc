// This is the version of the NLOX OLP interface which facilitates linking
// with FORTRAN programs. The function names have been modified, such that
// linking will succeed.
//
// To use the FORTRAN interface, you must add the line
//
//   include 'nlox_fortran_interface.f90'
//
// to any piece of FORTRAN code which will call these OLP functions.

#include <cstring>
#include <cmath>
#include "nlox_olp_fortran.h"

Process *fproc = NULL;

void nlox_olp_start_(char* fname, int* ierr) {
  fproc = new Process();
}

void nlox_olp_evalsubprocess_all_(char* sub, int* lsub, double* pp, int* next, double* mu, double* rval, double* acc) {
  size_t lsubs = *lsub
  std::string subs(sub,lsubs)
  int nextv = *next;
  double muv = *mu;
  fproc->evaluate_all(subs, pp, nextv,  muv, rval, acc);
}

void nlox_olp_evalsubprocess_(char* sub, int* lsub, char* typ, int* ltyp, char* cp, int* lcp, double* pp, int* next, double* mu, double* rval2, double* acc) {
  size_t lsubs = *lsub
  size_t ltyps = *ltyp;
  size_t lcps = *lcp;
  std::string subs(sub,lsubs)
  std::string typs(typ,ltyps);
  std::string cps(cp,lcps);
  int nextv = *next;
  double muv = *mu;
  if(cps.find("as")==0) fproc->evaluate_alpha(susb, typs, cps, pp, nextv, muv, rval2, acc);
  else fproc->evaluate(subs, typs, cps, pp, nextv, muv, rval2, acc);
}

void nlox_olp_evalsubprocess_cc_(char* sub, int* lsub, char* typ, int* ltyp, char* cp, int* lcp, double* pp, int* next, double* mu, double* rvalcc, double* acc){
  size_t lsubs = *lsub
  size_t ltyps = *ltyp;
  size_t lcps = *lcp;
  std::string subs(sub,lsubs)
  std::string typs(typ,ltyps);
  std::string cps(cp,lcps);
  int nextv = *next;
  double muv = *mu;
  
  if ( typs == "tree_loop" ) {
    std::cout << "Error:" << std::endl;
    std::cout << "Color correlated matrix elements are not generated for the specified intereference type \"" << typs << "\"." << std::endl;
    abort();
  }

  // For n=next external particles, the return array rvalcc will be of the form 
  // {Born, <T1T2>, <T1T3>,...,<T1Tn>,<T2T3>,<T2T4>,...,<T2Tn>,...,<T(n-1)Tn>}.

  std::string ccLStrings[7] = {"ccL0","ccL1","ccL2","ccL3","ccL4","ccL5","ccL6"};
  double rval2[3]; 

  // Set and evaluate the Born
  fproc->pc.set_param(ccLStrings[0],1.);
  for(int j=1; j<=nextv; j++) fproc->pc.set_param(ccLStrings[j],0.);
  if(cps.find("as")==0) fproc->evaluate_alpha(iv, typs, cps, pp, nextv, muv, rval2, acc);
  else  fproc->evaluate(iv, typs, cps, pp, nextv, muv, rval2, acc);
  rvalcc[0] = rval2[2];
  
  // Set and evaluate <TiTj>
  int cccounter = 1;
  fproc->pc.set_param(ccLStrings[0],0.);
  for(int j=1; j<=nextv; j++) {
    for(int k=j+1; k<=nextv; k++) {
      for(int l=0; l<=6; l++) {
        // Request <TjTk>=<TkTj>, k>j
        if(l==j||l==k) fproc->pc.set_param(ccLStrings[l],1.);
        else fproc->pc.set_param(ccLStrings[l],0.);
      }
      // Evaluate
      if(cps.find("as")==0) fproc->evaluate_alpha(iv, typs, cps, pp, nextv, muv, rval2, acc);
      else  fproc->evaluate(iv, typs, cps, pp, nextv, muv, rval2, acc);
      rvalcc[cccounter] = rval2[2];
      cccounter += 1;
    }
  }
}

void nlox_olp_setmass_(char* name, int* lname, double* mass, double* width, int* ierr ){
    if(fproc!=NULL){
        std::cout << "Warning: you are attempting to change a mass parameter after NLOX has been initialized,\n this change will not have an effect. Please set all mass parameters prior to the NLOX initialization call.\n";
        *ierr = 0;
    }
    else{
    size_t lnames = *lname;
    std::string mass_name(name,lnames);
    std::string width_name = mass_name;
    width_name.replace(0,1,"w");
    ProcessConst::set_param_in_file(mass_name,*mass,"nlox_parameters.par");
    ProcessConst::set_param_in_file(width_name,*width,"nlox_parameters.par");
    *ierr = 1;
    } 
}

void nlox_olp_setparameter_(char* para, int* lpar, double* re, double* im, int* ierr) {

  if (*im != 0.) { // No settable parameter is supported as complex in NLOX
    std::cout << "Warning: attempt to set unsupported imaginary part of parameter."
              << std::endl;
  }

  size_t lpars = *lpar;

  *ierr = fproc->pc.set_param(std::string(para,lpars), *re);
  fproc->update_mass(fproc->pc.mass_name(para));
}
