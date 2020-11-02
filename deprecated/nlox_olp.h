#ifndef __NLOX_OLP_H_
#define __NLOX_OLP_H_

#include "nlox_process.h"

#ifdef __cplusplus
extern "C" {
#endif

void NLOX_OLP_Start(char* fname, int* ierr);
void NLOX_OLP_EvalSubProcess_All(int* i, std::vector<FourVector> psp, int* next,double* mu, double* rval, double* acc);
void NLOX_OLP_EvalSubProcess(int* i, char* type, char* cp, std::vector<FourVector> psp, int* next, double* mu, double* rval2, double* acc);
void NLOX_OLP_EvalSubProcess_CC(int* i, char* type, char* cp, std::vector<FourVector> psp, int* next, double* Bij);
void NLOX_OLP_EvalSubProcess_SC(int* i, char* type, char* cp, std::vector<FourVector> psp, int* next, int pn, double** Bmunu);
void NLOX_OLP_EvalSubProcess_CC_and_SC(int* i, char* type, char* cp, std::vector<FourVector> psp, int* next, int pn, double*** Bmunuij);
void NLOX_OLP_SetMass(char* name, double* mass, double* width, int* ierr);
void NLOX_OLP_SetParameter(char* para, double* re, double* im, int* ierr);

#ifdef __cplusplus
}
#endif // __cplusplus

extern Process* proc;

#endif // __NLOX_OLP_H_
