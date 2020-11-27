#ifndef __KINEMATICS_H__
#define __KINEMATICS_H__

#include "Phase_Space_Tools.h"

namespace Kinematics{

    template<class T>
    T PseudoRapidity(FourVectorT<T> p){
        T pp = sqrt(p.p1*p.p1+p.p2*p.p2+p.p3*p.p3);
        T pz = p.p3;
        return 0.5*log((pp+pz)/(pp-pz));
    }
    
    template<class T>
    T Rapidity(FourVectorT<T> p){
        T E = p.p0;
        T pz = p.p3;
        return 0.5*log((E+pz)/(E-pz));
    }
    
    template<class T>
    T TransverseMomentum(FourVectorT<T> p){
        return sqrt(p.p1*p.p1+p.p2*p.p2);
    }
 
    template<class T>
    T RDistance(FourVectorT<T> p1, FourVectorT<T> p2){
        return 0;
    }

};

#endif
