#ifndef __KINEMATICS_H__
#define __KINEMATICS_H__

#include "Four_Vector.h"
#include <cmath>

namespace Kinematics{

    template<class T>
    T PseudoRapidity(FVectorT<T> p){
        T pp = sqrt(p.p1*p.p1+p.p2*p.p2+p.p3*p.p3);
        T pz = p.p3;
        return 0.5*log((pp+pz)/(pp-pz));
    }
    
    template<class T>
    T Rapidity(FVectorT<T> p){
        T E = p.p0;
        T pz = p.p3;
        return 0.5*log((E+pz)/(E-pz));
    }
    
    template<class T>
    T TransverseMomentum(FVectorT<T> p){
        return sqrt(p.p1*p.p1+p.p2*p.p2);
    }
 
    template<class T>
    T RDistance(FVectorT<T> p1, FVectorT<T> p2){
        T DEta = Rapidity(p1)-Rapidity(p2);
        T DPhi = std::atan(p1.p2/p1.p1) - std::atan(p2.p2/p2.p1);
        return std::sqrt(DEta*DEta + DPhi*DPhi);
    }

};

#endif
