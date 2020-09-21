#ifndef __NLOX_OLP_H_
#define __NLOX_OLP_H_

#include "nlox_process.h"

#ifdef __cplusplus
extern "C" {
#endif

void NLOX_OLP_Start(char* fname, int* ierr);
void NLOX_OLP_EvalSubProcess_All(char* SubProc, double* pp, int* next, double* mu, double* rval, double* acc);
void NLOX_OLP_EvalSubProcess(char* SubProc, char* type, char* cp,  double* pp, int* next, double* mu, double* rval2, double* acc);
void NLOX_OLP_EvalSubProcess_CC(char* SubProc, char* type, char* cp,  double* pp, int* next, double* mu, double* rvalcc, double* acc);
void NLOX_OLP_SetMass(char* name, double* mass, double* width, int* ierr);
void NLOX_OLP_SetParameter(char* para, double* re, double* im, int* ierr);

#ifdef __cplusplus
}
#endif // __cplusplus

extern Process* proc;

#endif // __NLOX_OLP_H_
