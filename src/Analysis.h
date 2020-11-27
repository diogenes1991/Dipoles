#ifndef __ANALYSIS_H__
#define __ANALYSIS_H__
#include "Phase_Space_Tools.h"
#include "Kinematics.h"

namespace Analysis{
    void ReweightEvent(FourVector* ExtMom, double* ExtMass, int* ExtPID, int BornNext, double* weight);
    void CombineParticles(FourVector* ExtMom, int* ExtPID, double ExtMass, int BorNext);
};

#endif
