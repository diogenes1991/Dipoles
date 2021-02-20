#ifndef __ANALYSIS_H__
#define __ANALYSIS_H__

#include "Kinematics.h"
#include "Utilities.h"

namespace Analysis{
    void ReweightEvent(FVector* ExtMom, double beta, double* ExtMass, int* ExtPID, int BornNext, double* weight);
    void InitializeHistograms(std::vector<Histogram>* Histograms);
    void FillHistograms(FVector* ExtMom, double beta, double* ExtMass, int* ExtPID, int BornNext, double xsec, std::vector<Histogram>* Histograms);
    void CombineParticles(FVector* ExtMom, double beta, int* ExtPID, double ExtMass, int BorNext);
};

#endif
