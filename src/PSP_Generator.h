#ifndef __GENERATOR_H_
#define __GENERATOR_H_

#include <vector>
#include <quadmath.h>
#include <cmath>
#include "Four_Vector.h"

#define ZERO_SPEED 1E-38

template <class T>
    T Lambda(T x, T y, T z){
    return x*x + y*y + z*z - 2*x*y - 2*y*z - 2*z*x;  
}

template <class T>
FMatrixT<T> Boost(T v_x, T v_y, T v_z){
  FMatrixT<T> out;
  T v = sqrt( v_x*v_x + v_y*v_y + v_z*v_z );
  FMatrixT<T> keta;
    if(v < ZERO_SPEED){
    keta.M[0][1] = 0 ;
    keta.M[0][2] = 0 ;
    keta.M[0][3] = 0 ;
    keta.M[1][0] = 0 ;
    keta.M[2][0] = 0 ;
    keta.M[3][0] = 0 ;
    }
      
    else{
    keta.M[0][1] = v_x / v ;
    keta.M[0][2] = v_y / v ;
    keta.M[0][3] = v_z / v ;
    keta.M[1][0] = v_x / v ;
    keta.M[2][0] = v_y / v ;
    keta.M[3][0] = v_z / v ;
    }
  
  T gamma = v*v;
    gamma = 1-gamma;
    gamma = sqrt(gamma);
    gamma = 1/gamma;
  
  out = 1 + gamma*v*keta + (gamma - 1)*keta*keta;
    
  return out;   
}

template <class T>
void Phase_Space_Point_Generator_2_Particle(FVectorT<T> P, FVectorT<T>& p1, T m1, FVectorT<T>& p2, T m2, T* r, T& J){
  
    T s = P.abs();
    T E = P.p0;

    T E1 = ( s + m1*m1 - m2*m2 ) / ( 2*sqrt(s) );
    T E2 = ( s + m2*m2 - m1*m1 ) / ( 2*sqrt(s) );
    T p  = sqrt(Lambda(s,m1*m1,m2*m2)/s)/2;  
    
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
    
    FMatrixT<T> B = Boost(vx,vy,vz);
    p1 = B * p1;
    p2 = B * p2;
    
    J = J*p;
    J = J/(4*M_PIq*sqrt(s));
}

template <class T>
void Recursive_PSP(FVectorT<T> P, int NPar, FVectorT<T>* PList, T* MList, T* r, T& J){
    
    //////////////////////////////////////////
    ////
    ////    PSP Generator based on:
    ////    Vernon and Barger: Collider Physics, 11.5
    ////    Its prescence here is to provide 
    ////    an arbitrary precision PSP
    ////    generator for testing.
    ////
    /////////////////////////////////////////

    FMatrixT<T> Boosts[NPar];
    for(int i=0;i<NPar;i++){Boosts[i] = 1;}
    
    int count = NPar;
    FVectorT<T> P_Left = P;
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
        
    FMatrixT<T> Boost_now;
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
    Phase_Space_Point_Generator_2_Particle(P_Left,PList[0],MList[0],PList[1],MList[1],ran,J);
    
    for(int i=0;i<NPar;i++){PList[i] = Boosts[i]*PList[i];}
}

template <class T>
void Repair_PSP(FVectorT<T>* PIn, T* MIn, int NPar, FVectorT<T>* POut, T* MOut){
        
    //////////////////////////////////////////////
    ////   
    //// This function corrects external moemnta to 
    //// the working precision of the class T.
    //// The argumets must be;
    ////   PIn: List with the two initial Four Vectors
    ////   MIn: List with the two initial masses 
    ////   NPar: the number os final particles
    ////   POut: List with the NPar final Four Vectors
    ////   MOut: List with the two initial masses 
    ////
    /////////////////////////////////////////////

    ////////////////////////////////////////
    ///
    ///  Boost into the CM of 
    ///  the initials to ensure the mappings 
    ///  are done correctly.
    ///
    ///////////////////////////////////////
    
    FVectorT<T> P = PIn[0]+PIn[1];        
    FMatrixT<T> BCM = Boost(-P.p1/P.p0,-P.p2/P.p0,-P.p3/P.p0);
    for(int i=0;i<2;i++)PIn[i]=BCM*PIn[i];
    for(int i=0; i<NPar;i++)POut[i] = BCM*POut[i]; 
    BCM = Boost(P.p1/P.p0,P.p2/P.p0,P.p3/P.p0);
    
    ////////////////////////////////////////
    ///
    ///   Finals-Square Correction 
    ///   Here we correct the value of m^2
    ///   for the final particles only
    ///
    /////////////////////////////////////////
    
    FVectorT<T> TP(0,0,0,0);
    for (int i=0;i<NPar;i++){
        FVectorT<T> PBar = POut[i].negateSpatial();
        T alpha = (MOut[i]*MOut[i] - POut[i]*POut[i]);
        T beta = PBar*POut[i];
        alpha = alpha/(beta + sqrt((beta*beta)+(POut[i]*POut[i]*alpha)));
        POut[i] = POut[i] + PBar*alpha;
        TP = TP - POut[i];
    }
    
    ////////////////////////////////////////
    ///
    ///   Final Sum Correction
    ///   After correcting the squares the 
    ///   total 3-momentum might not be zero 
    ///   hence we boost back into the frame 
    ///   where it is zero.
    ///
    /////////////////////////////////////////
    
    FMatrixT<T> B = Boost(-TP.p1/TP.p0,-TP.p2/TP.p0,-TP.p3/TP.p0);
    FVectorT<T> TPP(0,0,0,0);
    for (int i=0;i<NPar;i++){
        POut[i] = B*POut[i];
        TPP = TPP + POut[i];
    }
    
    ////////////////////////////////////////
    ///    
    ///   After correcting for the squares 
    ///   of the finals the total center of
    ///   mass energy might not match the 
    ///   initial center of mass energy.
    ///   We rescale the initials to this new Ecm.
    ///   Fisrt we pretend the initial masses 
    ///
    //////////////////////////////////////////////
    
    T c = 2*(PIn[0]*PIn[1]);
    T b = TPP*TPP - MIn[0]*MIn[0] - MIn[1]*MIn[1];
    T a = 0;
    FVectorT<T> PInB[2] = {PIn[0].negateSpatial(),PIn[1].negateSpatial()};
    T beta;
    
    b = b - (MIn[0]*MIn[0]*(PInB[0]*PIn[1]))/(PIn[0]*PInB[0]);
    b = b - (MIn[1]*MIn[1]*(PInB[1]*PIn[0]))/(PIn[1]*PInB[1]);
    b = b + ((MIn[1]*MIn[1]*(PIn[0]*PIn[0]))/(2*(PInB[0]*PIn[0])*(PInB[1]*PIn[1])))*(PInB[0]*PInB[1]);
    b = b + ((MIn[0]*MIn[0]*(PIn[1]*PIn[1]))/(2*(PInB[1]*PIn[1])*(PInB[0]*PIn[0])))*(PInB[1]*PInB[0]);
    b = b + ((MIn[0]*MIn[0])*(PIn[0]*PIn[0])*(PIn[0]*PIn[0]))/(2*(PIn[0]*PInB[0])*(PIn[0]*PInB[0]));
    b = b + ((MIn[1]*MIn[1])*(PIn[1]*PIn[1])*(PIn[1]*PIn[1]))/(2*(PIn[1]*PInB[1])*(PIn[1]*PInB[1]));
    
    c = c - ((PIn[0]*PIn[0])*(PInB[0]*PIn[1]))/(PIn[0]*PInB[0]);
    c = c - ((PIn[1]*PIn[1])*(PInB[1]*PIn[0]))/(PIn[1]*PInB[1]);
    c = c + ((PInB[0]*PInB[1])*(PIn[0]*PIn[0])*(PIn[1]*PIn[1]))/(2*(PIn[0]*PInB[0])*(PIn[1]*PInB[1]));
    c = c + ((PIn[0]*PIn[0])*(PIn[0]*PIn[0])*(PIn[0]*PIn[0]))/(4*(PIn[0]*PInB[0])*(PIn[0]*PInB[0]));
    c = c + ((PIn[1]*PIn[1])*(PIn[1]*PIn[1])*(PIn[1]*PIn[1]))/(4*(PIn[1]*PInB[1])*(PIn[1]*PInB[1]));
    
    a = a + (MIn[0]*MIn[0]*MIn[0]*MIn[0]*(PIn[0]*PIn[0]))/(4*(PIn[0]*PInB[0])*(PIn[0]*PInB[0]));
    a = a + (MIn[1]*MIn[1]*MIn[1]*MIn[1]*(PIn[1]*PIn[1]))/(4*(PIn[1]*PInB[1])*(PIn[1]*PInB[1]));
    a = a + (MIn[0]*MIn[0]*MIn[1]*MIn[1]*(PInB[0]*PInB[1]))/(2*(PIn[0]*PInB[0])*(PIn[1]*PInB[1]));
    
    beta = (2*c)/(b+sqrt((b*b)-(4*a*c)));
            
    for(int i=0;i<2;i++){
        FVectorT<T> PBar= PIn[i].negateSpatial();
        T alpha = (MIn[i]*MIn[i]*beta - PIn[i]*PIn[i]);
        alpha = alpha/(2*(PBar*PIn[i]));
        PIn[i] = PIn[i] + PBar*alpha;
        PIn[i] = PIn[i]*(1/sqrt(beta));
    }
    
    P = PIn[0]+PIn[1];
    B = Boost(-P.p1/P.p0,-P.p2/P.p0,-P.p3/P.p0);
    PIn[0] = B*PIn[0];
    PIn[1] = B*PIn[1];
    FVectorT<T> PP=B*P;
    
    for(int i=0;i<2;i++)PIn[i]=BCM*PIn[i];
    for(int i=0; i<NPar;i++)POut[i] = BCM*POut[i];
}

#endif

