// This is the version of the NLOX OLP interface which facilitates linking
// with FORTRAN programs. The function names have been modified, such that
// linking will succeed.
//
// To use the FORTRAN interface, you must add the line
//
//   include 'nlox_fortran_interface.f90'
//
// to any piece of FORTRAN code which will call these OLP functions.

#ifndef __NLOX_OLP_FORTRAN_H_
#define __NLOX_OLP_FORTRAN_H_

#include "nlox_process.h"

#ifdef __cplusplus
extern "C" {
#endif

void nlox_olp_start_(char* fname, int* ierr);
void nlox_olp_evalsubprocess_all_(int* i, double* pp, int* next, double* mu, double* rval, double* acc);
void nlox_olp_evalsubprocess_(int* i, char* typ, int* ltyp, char* cp, int* lcp, double* pp, int* next, double* mu, double* rval2, double* acc);
void nlox_olp_evalsubprocess_cc_(int* i, char* typ, int* ltyp, char* cp, int* lcp, double* pp, int* next, double* mu, double* rvalcc, double* acc);
void nlox_olp_setmass_(char* name, int* lname, double* mass, double* width, int* ierr );
void nlox_olp_setparameter_(char* para, int* lpar, double* re, double* im, int* ierr);

#ifdef __cplusplus
}
#endif // __cplusplus

extern Process* fproc;

#endif // __NLOX_OLP_FORTRAN_H_
