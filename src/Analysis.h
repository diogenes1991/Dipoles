#ifndef __ANALYSIS_H__
#define __ANALYSIS_H__

#include "Kinematics.h"

namespace Analysis{
    void ReweightEvent(FVector* ExtMom, double* ExtMass, int* ExtPID, int BornNext, double* weight);
    void CombineParticles(FVector* ExtMom, int* ExtPID, double ExtMass, int BorNext);
};

#endif
