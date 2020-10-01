#ifndef __GENERATOR_H_
#define __GENERATOR_H_

#include <vector>
#include "Phase_Space_Tools.h"

template <class T>
void Base_Case(FourVectorT<T> P, FourVectorT<T>& p1, T m1, FourVectorT<T>& p2, T m2, T* r, T& J){
  
    T s = P.abs();
    T E = P.p0;

    T E1 = ( s + m1*m1 - m2*m2 ) / ( 2*sqrt(s) );
    T E2 = ( s + m2*m2 - m1*m1 ) / ( 2*sqrt(s) );
    T p  = sqrt(lambda(s,m1*m1,m2*m2)/s)/2;  
    
    T cos_th = 2 * r[0] - 1 ;
    T ph = 2*M_PIq*r[1];
    T sin_th = sqrt(1 - cos_th*cos_th);
    
    p1.p0 =  E1;
    p1.p1 =  -p*cos(ph)*sin_th;
    p1.p2 =  -p*sin(ph)*sin_th;
    p1.p3 =  -p*cos_th;
    
    p2.p0 =  E2;
    p2.p1 =  p*cos(ph)*sin_th;
    p2.p2 =  p*sin(ph)*sin_th;
    p2.p3 =  p*cos_th;
    
    T vx = P.p1 / E;
    T vy = P.p2 / E;
    T vz = P.p3 / E;
    
    FourMatrixT<T> B = Boost(vx,vy,vz);
    p1 = B * p1;
    p2 = B * p2;
    
    J = J*p;
    J = J/(4*M_PIq*sqrt(s));
}

template <class T>
void Recursive_PSP_Generator(FourVectorT<T> P, std::vector<FourVectorT<T>> PList, std::vector<T> MList, T* r, T& J){
    
    //////////////////////////////////////////
    ////
    ////    PSP Generator based on:
    ////    Vernon and Barger: Collider Physics, 11.5
    ////    Its prescence here is to provide 
    ////    an arbitrary precision PSP
    ////    generator for testing.
    ////
    /////////////////////////////////////////

    int NPar = PList.size();

    FourMatrixT<T> Boosts[NPar];
    for(int i=0;i<NPar;i++){Boosts[i] = 1;}
    
    int count = NPar;
    FourVectorT<T> P_Left = P;
    while (count > 2){
    T s = P_Left.abs();
    T E = P_Left.p0;
    
    T s_min = 0;
    for(int i=0;i<(count-1);i++){ s_min += MList[i]; }
    s_min = s_min*s_min;
    T s_max = (sqrt(s) - MList[count-1])*(sqrt(s) - MList[count-1]);
    
    T s_now = (s_max-s_min)*r[3*(NPar-count)] + s_min;
    
    T E_now = (s - s_now + MList[count-1]*MList[count-1]) / (2*sqrt(s)) ;
    T pp_now = sqrt( E_now*E_now - MList[count-1]*MList[count-1] );
    T cos_th_now = 2 * r[3*(NPar-count)+1] - 1  ;
    T sin_th_now = sqrt( 1 - cos_th_now*cos_th_now );
    T ph_now = 2*M_PIq*r[3*(NPar-count)+2] ;
    
    J = J*pp_now;
    J = J*(s_max-s_min);
    J = J/(8*M_PIq*M_PIq*sqrt(s));
    
    PList[count-1].p0 =  E_now ;
    PList[count-1].p1 =  pp_now*cos(ph_now)*sin_th_now ;
    PList[count-1].p2 =  pp_now*sin(ph_now)*sin_th_now ;
    PList[count-1].p3 =  pp_now*cos_th_now ;
        
    FourMatrixT<T> Boost_now;
    Boost_now = 1;
    
    
    T vx = P_Left.p1 / E;
    T vy = P_Left.p2 / E;
    T vz = P_Left.p3 / E;  
    Boost_now = Boost(vx,vy,vz);
    
    for(int i=0;i<=(count-1);i++)Boosts[i] = Boosts[i]*Boost_now;
    
    P_Left.p0=sqrt(s)-PList[count-1].p0;
    P_Left.p1=-PList[count-1].p1;
    P_Left.p2=-PList[count-1].p2;
    P_Left.p3=-PList[count-1].p3;
    
    count = count - 1;

    }
    
    T ran[2] = {r[3*NPar-6],r[3*NPar-5]};
    Base_Case(P_Left,PList[0],MList[0],PList[1],MList[1],ran,J);
    
    for(int i=0;i<NPar;i++){PList[i] = Boosts[i]*PList[i];}
}

#define Get_PSP Recursive_PSP_Generator

#endif

