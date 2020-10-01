#ifndef _Dipoles_Definitions_H
#define _Dipoles_Definitions_H

#include "Phase_Space_Tools.h"
// #include "Plus_Distribution.h"

#include <gsl/gsl_sf_dilog.h>
#define PolyLog2 gsl_sf_dilog

// 
//  The log_4pi definition is needed if the NLO Cross-Sections
//  are defined not to have (4pi)^epsilon/Gamma(1-epsilon) upfront
//  Otherwise this log is never used.
//

#define log_4pi 2.5310242469692907929778915942694118477
#define pi2 9.869604401089358618834490999876151135

// 
//   II Functions and Maps 
// 

template <class T>
T x_ab(FourVectorT<T> pa, FourVectorT<T> pb, FourVectorT<T> k){
    T out = pa*pb - pa*k - pb*k;
    out = out / (pa*pb);
    return out;
}

template <class T>
T xab_cut(FourVectorT<T> pa, FourVectorT<T> pb, T ma, T mb){ 
    return ((ma*mb)/(pa*pb));
}

template <class T>
T Rab(FourVectorT<T> pa, FourVectorT<T> pb, T x){
    T aux = pow(2*(pa*pb)*x,2)-((pa*pa)*(pb*pb));
    aux = aux / (pow(2*(pa*pb),2)-((pa*pa)*(pb*pb)));
    return sqrt(aux);
}

template <class T>
FourVectorT<T> pa_II(FourVectorT<T> pa, FourVectorT<T> pb, T xab){
    T sbar = 2*(pa*pb);
    FourVectorT<T> aux = (2*sbar*(pa*pa)*(1.0-xab*xab))*pb;
    aux = 1.0/((Rab(pa,pb,xab)+xab)*(sbar*sbar-4*(pa*pa)*(pb*pb)))*aux;
    aux = Rab(pa,pb,xab)*pa + aux;
    return aux;
}

template <class T>
FourVectorT<T> pa_II(FourVectorT<T> pa, FourVectorT<T> pb, FourVectorT<T> k){
    return pa_II(pa,pb,x_ab(pa,pb,k));
}

template <class T>
FourVectorT<T> Boost_II(FourVectorT<T> pa, FourVectorT<T> pb, FourVectorT<T> k, FourVectorT<T> ki){
    FourVectorT<T> Pab = pa + pb - k;
    FourVectorT<T> Pab_tilde = pa_II(pa,pb,k) + pb;
    FourVectorT<T> out = Pab_tilde;
    out = out * (2*(Pab*ki)/(Pab*Pab));
    FourVectorT<T> aux = Pab+Pab_tilde;
    aux = ((aux*ki)/(aux*Pab))*aux;
    return ki - aux + out;
}

template <class T>
void Build_II_Momenta(int NEXT_BORN, std::vector<FourVectorT<T>> P_DAT, std::vector<FourVectorT<T>> P_TIL, int EMIT, int SPEC, int RADI){
    P_TIL.at(EMIT) = pa_II(P_DAT.at(EMIT),P_DAT.at(SPEC),P_DAT.at(RADI));
    for(int II=2;II<NEXT_BORN;II++){
        P_TIL.at(II) = Boost_II(P_DAT.at(EMIT),P_DAT.at(SPEC),P_DAT.at(RADI),P_DAT.at(II));
    }
}

template <class T>
T g_ab_fermion(FourVectorT<T> pa, FourVectorT<T> pb, FourVectorT<T> k, T ma, T mb){
    T xab = x_ab(pa,pb,k);
    if (xab<xab_cut(pa,pb,ma,mb)) return 0;
    T out = 2.0/(1.0-xab);
    out  = out - (1.0+xab);
    out  = out - (xab*(pa*pa)/(pa*k));
    out  = out / (xab*(pa*k));
    return out;   
}

// NOTE: Currently g_ab_boson is a placeholder and its overloaded to g_ab_fermion temporarily //
template <class T>
T g_ab_boson(FourVectorT<T> pa, FourVectorT<T> pb, FourVectorT<T> k){
    
    return g_ab_fermion(pa,pb,k);    
}

// template <class T>
// T CurlyG_ab_fermion(FourVectorT<T> pa, FourVectorT<T> pb, T ma, T mb, T x, T Ix, T I1){
    
//     T sab = ma*ma + mb*mb + 2*(pa*pb);
//     T sab_bar = 2*(pa*pb);
//     T LOG_DR = /*log_4pi + */log(mu*mu*sab/(sab_bar*sab_bar));



// };

template <class T>
void G_ab_fermion(FourVectorT<T> pa, FourVectorT<T> pb, T ma, T mb, T mu, T* RVAL){

    // TODO: This template only works for T=double due to the gsl_sf_dilog call
 
    T sab = ma*ma + mb*mb + 2*(pa*pb);
    T sab_bar = 2*(pa*pb);
    T x0 = xab_cut(pa,pb,ma,mb);
    T LOG_DR = /*log_4pi + */log(mu*mu*sab/(sab_bar*sab_bar)) - 2*log(1.0-x0);

    if (ma==0.0){

        RVAL[0] = 1;
        RVAL[1] = - 3.0/2.0 - LOG_DR - x0*(1.0+x0)/2.0;
        RVAL[2] = LOG_DR*LOG_DR - pi2 / 4;
    }

    else{

        T BETA = sab_bar + 2*ma*ma - sqrt(lambda(sab,ma*ma,mb*mb));
          BETA = BETA / (sab_bar + 2*ma*ma + sqrt(lambda(sab,ma*ma,mb*mb)));
        T  LOG = log(BETA);
        
        RVAL[0] = 0;
        RVAL[1] = - 1 - sab_bar/lambda(sab,ma*ma,mb*mb)*LOG;
        RVAL[2] = (2*PolyLog2(1-BETA) + (1.0/2.0)*LOG*LOG + ((sab_bar+2*ma*ma)/(sab_bar))*LOG);
        RVAL[2] = LOG_DR*RVAL[1]-(sab_bar/lambda(sab,ma*ma,mb*mb))*RVAL[2];

    }

}

//
//  IF & FI Functions and Maps 
//

template <class T>
T x_ia(FourVectorT<T> pa, FourVectorT<T> pi, FourVectorT<T> k){
    
    return ((pa*pi) + (pa*k) - (pi*k)) / ( (pa*pi) + (pa*k) );
}

template <class T>
T z_ia(FourVectorT<T> pa, FourVectorT<T> pi, FourVectorT<T> k){
    
    return (pa*pi) / ( (pa*pi) + (pa*k) );
}

template <class T>
T R_ia(FourVectorT<T> pa, FourVectorT<T> pi, FourVectorT<T> k, T x){
    FourVectorT<T> P_ai = pi - pa + k;
    T P_ai_bar = (P_ai*P_ai) - (pa*pa) - (pi*pi);
    
    T out = (P_ai_bar + 2*(pa*pa)*x)*(P_ai_bar + 2*(pa*pa)*x);
    out = out - 4*(pa*pa)*(P_ai*P_ai)*x*x;
    out = out / (lambda((P_ai*P_ai),(pa*pa),(pi*pi)));
    out = sqrt(out);
    return out;  
}

template <class T>
T g_ai_fermion(FourVectorT<T> pa, FourVectorT<T> pi, FourVectorT<T> k){
  
    T xai = x_ia(pa,pi,k);
    T zai = z_ia(pa,pi,k);
  
    T Ria = R_ia(pa,pi,k,xai);   
  
    T aux1 = 2.0 / (2.0 - xai - zai );
    T aux2 = (1.0 + xai) * Ria;
    T aux3 = xai*(pa*pa)/(k*pa);
  
    T out = aux1 - aux2 - aux3;
    out = out / ((pa*k)*xai);
    return out;
}

template <class T>
T g_ia_fermion(FourVectorT<T> pi, FourVectorT<T> pa, FourVectorT<T> k){
    
    T xia = x_ia(pa,pi,k);
    T zia = z_ia(pa,pi,k);
    
    T out = 2./(2.-xia-zia);
    out = out - 1 - zia - (pi*pi)/(pi*k);
    out = out/(xia*(pi*k));
    
    return out;   
}

// NOTE: Currently g_ia(ai)_boson are a placeholders and are overloaded to g_ia(ai)_fermion temporarily //
template <class T>
T g_ai_boson(FourVectorT<T> pa, FourVectorT<T> pi, FourVectorT<T> k){
    
    return g_ai_fermion(pa,pi,k);    
}

template <class T>
T g_ia_boson(FourVectorT<T> pi, FourVectorT<T> pa, FourVectorT<T> k){
    
    return g_ia_fermion(pi,pa,k);
}

template <class T>
FourVectorT<T> pi_IFFI(FourVectorT<T> pa, FourVectorT<T> pi, FourVectorT<T> k){
  FourVectorT<T> P_ia = pi + k - pa;
  FourVectorT<T> out1 = pa - ( P_ia*pa / ( P_ia * P_ia ) ) * P_ia ;
  FourVectorT<T> out2 = ( ( P_ia*P_ia - pa*pa + pi*pi ) / (2*(P_ia*P_ia) ) ) * P_ia ;
  FourVectorT<T> aux = (pi+k);
  T prefactor = sqrt( lambda(P_ia*P_ia,pa*pa,pi*pi) ) / sqrt( lambda(aux*aux,P_ia*P_ia,pa*pa));
  return prefactor*out1 + out2;  
}

template <class T>
FourVectorT<T> pa_IFFI(FourVectorT<T> pa, FourVectorT<T> pi, FourVectorT<T> k){
  
  return pi_IFFI(pa,pi,k)- pi - k + pa ;  
}

template <class T>
void Build_IF_Momenta(int NEXT_BORN, std::vector<FourVectorT<T>> P_DAT, std::vector<FourVectorT<T>> P_TIL, int EMIT, int SPEC, int RADI){
    P_TIL.at(EMIT) = pa_IFFI(P_DAT.at(EMIT),P_DAT.at(SPEC),P_DAT.at(RADI));
    P_TIL.at(SPEC) = pi_IFFI(P_DAT.at(EMIT),P_DAT.at(SPEC),P_DAT.at(RADI));
}

template <class T>
void Build_FI_Momenta(int NEXT_BORN, std::vector<FourVectorT<T>> P_DAT, std::vector<FourVectorT<T>> P_TIL, int EMIT, int SPEC, int RADI){
    P_TIL.at(EMIT) = pi_IFFI(P_DAT.at(EMIT),P_DAT.at(SPEC),P_DAT.at(RADI));
    P_TIL.at(SPEC) = pa_IFFI(P_DAT.at(EMIT),P_DAT.at(SPEC),P_DAT.at(RADI));
}

// 
//  FF Functions and Maps
//

template <class T>
T y_ij(FourVectorT<T> pi,FourVectorT<T> pj, FourVectorT<T> k){
    T aux = pi*k;
           aux = aux / (pi*pj+pi*k+pj*k);
    return aux;
}

template <class T>
T z_ij(FourVectorT<T> pi,FourVectorT<T> pj, FourVectorT<T> k){
    T aux = pi*pj;
           aux = aux / (pi*pj+pj*k);
    return aux;
}

template <class T>
T R_ij(FourVectorT<T> pi, FourVectorT<T> pj, FourVectorT<T> k,T y){
    FourVectorT<T> Pij = pi + pj + k;
    T Pijbar = Pij*Pij - pi*pi - pj*pj;
    T aux = (2*(pj*pj)-Pijbar*(1.0-y));
           aux = aux*aux;
           aux = aux - 4*(Pij*Pij)*(pj*pj);
           aux = aux / lambda(Pij*Pij,pi*pi,pj*pj);
           aux = sqrt(aux);
           return aux;
}

template <class T>
T g_ij_fermion(FourVectorT<T> pi, FourVectorT<T> pj, FourVectorT<T> k){
    T zij = z_ij(pi,pj,k);
    T yij = y_ij(pi,pj,k);
    T aux = zij*(1.0-yij);
           aux = 1.0 - aux;
           aux = 2.0 / aux;
           aux = aux - (1.0+zij);
           aux = aux - (pi*pi)/(pi*k);
           aux = aux / (R_ij(pi,pj,k,yij)*(pi*k));
           return aux;      
}

template <class T>
T g_ij_boson(FourVectorT<T> pi, FourVectorT<T> pj, FourVectorT<T> k){
    return g_ij_fermion(pi,pj,k);  
}

template <class T>
FourVectorT<T> pj_FF(FourVectorT<T> pi, FourVectorT<T> pj, FourVectorT<T> k){
    FourVectorT<T> Pij = pi + pj + k;
    T yij = y_ij(pi,pj,k);
    FourVectorT<T> aux = ((Pij*pj)/(Pij*Pij))*Pij;
           aux = pj - aux;
           aux = aux *(1.0/R_ij(pi,pj,k,yij));
           aux = aux + ((Pij*Pij+pj*pj-pi*pi)/(2*(Pij*Pij)))*Pij;
           return aux;           
}

template <class T>
FourVectorT<T> pi_FF(FourVectorT<T> pi, FourVectorT<T> pj, FourVectorT<T> k){
    return pi + pj + k - pj_FF(pi,pj,k);
}

template <class T>
void Build_FF_Momenta(int NEXT_BORN, std::vector<FourVectorT<T>> P_DAT, std::vector<FourVectorT<T>> P_TIL, int EMIT, int SPEC, int RADI){
    P_TIL.at(EMIT) = pi_FF(P_DAT.at(EMIT),P_DAT.at(SPEC),P_DAT.at(RADI));
    P_TIL.at(SPEC) = pj_FF(P_DAT.at(EMIT),P_DAT.at(SPEC),P_DAT.at(RADI));
    for(int II=0;II<NEXT_BORN;II++){
        if ( II!=EMIT && II!=SPEC ){
            P_TIL.at(II) = P_DAT.at(II);
        }
    }
}

#endif
