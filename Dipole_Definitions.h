#ifndef _Dipoles_Definitions_H
#define _Dipoles_Definitions_H

#include "Integrand_Classes.h"


// II Functions and Maps ///

double x_ab(Vector pa, Vector pb, Vector k){
    double out = pa*k;
           out = out + pb*k;
           out = out / (pa*pb);
           out = 1.-out;
           return out;
}

double Rab(Vector pa, Vector pb, double x){
    double aux = pow(2,(pa*pb)*x)-((pa*pa)*(pb*pb));
           aux = aux / (pow(2,(pa*pb))-((pa*pa)*(pb*pb)));
           return sqrt(aux);
}

double Rab(double s, double ma, double mb, double x){
    return sqrt(((s-ma*ma-mb*mb)*(s-ma*ma-mb*mb)*x*x-4*ma*ma*mb*mb)/lambda(s,ma*ma,mb*mb));
    
}

Vector pa_II(Vector pa, Vector pb, Vector k){
    double sbar = 2*(pa*pb);
    Vector aux = pa - (sbar/(2*(pb*pb)))*pb;
           aux = Rab(pa,pb,x_ab(pa,pb,k))*aux;
           aux = aux - (sbar*x_ab(pa,pb,k)/(2*(pb*pb)))*pb;
           return aux;
}

Vector Boost_II(Vector pa, Vector pb, Vector k, Vector ki){
    Vector Pab = pa + pb -k;
    Vector Pab_tilde = pa_II(pa,pb,k) + pb;
    Vector out = Pab_tilde;
           out = out * (2*(Pab*ki)/(Pab*Pab));
    Vector aux = Pab+Pab_tilde;
           aux = ((Pab*ki+Pab_tilde*ki)/(Pab*Pab+Pab_tilde*Pab))*aux;
           return out + aux;

}

void Build_II_Momenta(int NEXT_BORN, Vector* P_DAT, Vector* P_TIL, int EMIT, int SPEC){
    P_TIL[EMIT-1] = pa_II(P_DAT[EMIT-1],P_DAT[SPEC-1],P_DAT[NEXT_BORN]);
    P_TIL[SPEC-1] = P_DAT[SPEC-1];
    for(int II=2;II<NEXT_BORN;II++){
        P_TIL[II] = Boost_II(P_DAT[EMIT-1],P_DAT[SPEC-1],P_DAT[NEXT_BORN],P_DAT[II]);
    }
    
    
}

double g_ab_fermion(Vector pa, Vector pb, Vector k){
    
    double xab = x_ab(pa,pb,k);
    double out = 2./(1-xab);
           out  = out - (1.+xab);
           out  = out - (xab*(pa*pa)/(pa*k));
           out  = out / (xab*(pa*k));
    return out;
    
}

// NOTE: Currently g_ab_boson is a placeholder and its overloaded to g_ab_fermion temporarily //

double g_ab_boson(Vector pa, Vector pb, Vector k){
    
    return g_ab_fermion(pa,pb,k);
    
}





Plus_Distribution CurlyG_ab_fermion();






// IF & FI Functions and Maps //


double x_ia(Vector pa, Vector pi, Vector k){
    return ((pa*pi) + (pa*k) - (pi*k)) / ( (pa*pi) + (pa*k) );
}

double z_ia(Vector pa, Vector pi, Vector k){
    return (pa*pi) / ( (pa*pi) + (pa*k) );
}

double R_ia(Vector pa, Vector pi, Vector k, double x){
    Vector P_ai = pi - pa + k;
    double P_ai_bar = (P_ai*P_ai) - (pa*pa) - (pi*pi);
    
    double out = (P_ai_bar + 2*(pa*pa)*x)*(P_ai_bar + 2*(pa*pa)*x);
    out = out - 4*(pa*pa)*(P_ai*P_ai)*x*x;
    out = out / (lambda((P_ai*P_ai),(pa*pa),(pi*pi)));
    out = sqrt(out);
    return out;
    
}

double g_ai_fermion(Vector pa, Vector pi, Vector k){
  
    double xai = x_ia(pa,pi,k);
    double zai = z_ia(pa,pi,k);
  
    double Ria = R_ia(pa,pi,k,xai);   
  
    double aux1 = 2.0 / (2.0 - xai - zai );
    double aux2 = (1.0 + xai) * Ria;
    double aux3 = xai*(pa*pa)/(k*pa);
  
    double out = aux1 - aux2 - aux3;
    out = out / ((pa*k)*xai);
    return out;
  
  
}

double g_ia_fermion(Vector pi, Vector pa, Vector k){
    
    double xia = x_ia(pa,pi,k);
    double zia = z_ia(pa,pi,k);
    
    double out = 2./(2.-xia-zia);
    out = out - 1 - zia - (pi*pi)/(pi*k);
    out = out/(xia*(pi*k));
    
    return out;
    
    
    
}

// NOTE: Currently g_ia(ai)_boson are a placeholders and are overloaded to g_ia(ai)_fermion temporarily //

double g_ai_boson(Vector pa, Vector pi, Vector k){
    
    return g_ai_fermion(pa,pi,k);
    
}

double g_ia_boson(Vector pi, Vector pa, Vector k){
    return g_ia_fermion(pi,pa,k);
}

Vector pi_IFFI(Vector pa, Vector pi, Vector k){
  Vector P_ia = pi + k - pa;
  Vector out1 = pa - ( P_ia*pa / ( P_ia * P_ia ) ) * P_ia ;
  Vector out2 = ( ( P_ia*P_ia - pa*pa + pi*pi ) / (2*(P_ia*P_ia) ) ) * P_ia ;
  Vector aux = (pi+k);
  double prefactor = sqrt( lambda(P_ia*P_ia,pa*pa,pi*pi) ) / sqrt( lambda(aux*aux,P_ia*P_ia,pa*pa));
  return prefactor*out1 + out2;  
}

Vector pa_IFFI(Vector pa, Vector pi, Vector k){
  return pi_IFFI(pa,pi,k)- pi - k + pa ;  
}

void Build_IF_Momenta(int NEXT_BORN, Vector* P_DAT, Vector* P_TIL, int EMIT, int SPEC){
    P_TIL[EMIT-1] = pa_IFFI(P_DAT[EMIT-1],P_DAT[SPEC-1],P_DAT[NEXT_BORN]);
    P_TIL[SPEC-1] = pi_IFFI(P_DAT[EMIT-1],P_DAT[SPEC-1],P_DAT[NEXT_BORN]);
    for(int II=0;II<NEXT_BORN;II++){
        P_TIL[II] = P_DAT[II];
    }
}

void Build_FI_Momenta(int NEXT_BORN, Vector* P_DAT, Vector* P_TIL, int EMIT, int SPEC){
    P_TIL[EMIT-1] = pi_IFFI(P_DAT[EMIT-1],P_DAT[SPEC-1],P_DAT[NEXT_BORN]);
    P_TIL[SPEC-1] = pa_IFFI(P_DAT[EMIT-1],P_DAT[SPEC-1],P_DAT[NEXT_BORN]);
    for(int II=0;II<NEXT_BORN;II++){
        if ( II!=(EMIT-1) && II!=(SPEC-1) ){
            P_TIL[II] = P_DAT[II];
        }
    }
}


// FF Functions and Maps //

double y_ij(Vector pi,Vector pj, Vector k){
    double aux = pi*k;
           aux = aux / (pi*pj+pi*k+pj*k);
    return aux;
}

double z_ij(Vector pi,Vector pj, Vector k){
    double aux = pi*pj;
           aux = aux / (pi*pj+pj*k);
    return aux;
}

double R_ij(Vector pi, Vector pj, Vector k,double y){
    Vector Pij = pi + pj + k;
    double Pijbar = Pij*Pij - pi*pi - pj*pj;
    double aux = (2*(pj*pj)-Pijbar*(1.0-y));
           aux = aux*aux;
           aux = aux - 4*(Pij*Pij)*(pj*pj);
           aux = aux / lambda(Pij*Pij,pi*pi,pj*pj);
           aux = sqrt(aux);
           return aux;
}

double g_ij_fermion(Vector pi, Vector pj, Vector k){
    double zij = z_ij(pi,pj,k);
    double yij = y_ij(pi,pj,k);
    double aux = zij*(1.0-yij);
           aux = 1.0 - aux;
           aux = 2.0 / aux;
           aux = aux - (1.0+zij);
           aux = aux - (pi*pi)/(pi*k);
           aux = aux / (R_ij(pi,pj,k,yij)*(pi*k));
           return aux;
    
    
}

double g_ij_boson(Vector pi, Vector pj, Vector k){
    
    return g_ij_fermion(pi,pj,k);
    
}


Vector pj_FF(Vector pi, Vector pj, Vector k){
    Vector Pij = pi + pj + k;
    double yij = y_ij(pi,pj,k);
    Vector aux = ((Pij*pj)/(Pij*Pij))*Pij;
           aux = pj - aux;
           aux = aux *(1.0/R_ij(pi,pj,k,yij));
           aux = aux + ((Pij*Pij+pj*pj-pi*pi)/(2*(Pij*Pij)))*Pij;
           return aux;
             
}

Vector pi_FF(Vector pi, Vector pj, Vector k){
    return pi + pj + k - pj_FF(pi,pj,k);
}

void Build_FF_Momenta(int NEXT_BORN, Vector* P_DAT, Vector* P_TIL, int EMIT, int SPEC){
    P_TIL[EMIT-1] = pi_FF(P_DAT[EMIT-1],P_DAT[SPEC-1],P_DAT[NEXT_BORN]);
    P_TIL[SPEC-1] = pj_FF(P_DAT[EMIT-1],P_DAT[SPEC-1],P_DAT[NEXT_BORN]);
    for(int II=0;II<NEXT_BORN;II++){
        if ( II!=(EMIT-1) && II!=(SPEC-1) ){
            P_TIL[II] = P_DAT[II];
        }
    }
}



#endif
