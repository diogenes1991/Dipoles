#ifndef _Phase_Spaces_Real_H
#define _Phase_Spaces_Real_H
#include <stdlib.h>

const double PI = 3.1415926535897932385;

double aleatorio(double a, double b){
  double r = float(rand())/RAND_MAX;
  
  return a + (b-a)*r;
  
  
  
}

double lambda(double x, double y, double z){
  return x*x + y*y + z*z - 2*x*y - 2*y*z - 2*z*x;  
}

CMatrix Boost(double v_x, double v_y, double v_z){
  CMatrix out;
  
  double v = sqrt( v_x*v_x + v_y*v_y + v_z*v_z );
  
//   cout << "Boost Triggered with v = ( " << v_x << " , " << v_y << " , " << v_z << " )" << endl;  
  
  CMatrix keta;
    if(v < 1E-9)
    {
    keta.M[0][1] = 0.0 ;
    keta.M[0][2] = 0.0 ;
    keta.M[0][3] = 0.0 ;
    keta.M[1][0] = 0.0 ;
    keta.M[2][0] = 0.0 ;
    keta.M[3][0] = 0.0 ;
    }
      
    else
    {
    keta.M[0][1] = v_x / v ;
    keta.M[0][2] = v_y / v ;
    keta.M[0][3] = v_z / v ;
    keta.M[1][0] = v_x / v ;
    keta.M[2][0] = v_y / v ;
    keta.M[3][0] = v_z / v ;
    }
  
  double gamma = sqrt(1.0/(1.0 - v*v));
  
  out = 1.0 + gamma*v*keta + (gamma - 1.0)*keta*keta;
  
//   cout << "Boost subroutine triggered with \nB = \n" << out << endl;
  
  return out;
    
}

void Phase_Space_Point_Generator_2_Particle(Vector P, Vector& p1, double m1, Vector& p2, double m2, double* r, double& J){
  
  double s = P*P;
  double E = P.get_t();
  
  double E1 = ( s + m1*m1 - m2*m2 ) / ( 2*sqrt(s) );
  double E2 = ( s + m2*m2 - m1*m1 ) / ( 2*sqrt(s) );
  double p  = (1.0/2.0)*sqrt(lambda(s,m1*m1,m2*m2)/s);  
  
  double cos_th = 2 * r[0] - 1.0 ;
  double ph = 2.0*PI*r[1];
  double sin_th = sqrt(1.0 - cos_th*cos_th);
  
  p1.set_t( E1 );
  p1.set_x( -p*cos(ph)*sin_th );
  p1.set_y( -p*sin(ph)*sin_th );
  p1.set_z( -p*cos_th );
  
  p2.set_t( E2 );
  p2.set_x( p*cos(ph)*sin_th );
  p2.set_y( p*sin(ph)*sin_th );
  p2.set_z( p*cos_th );
  
  ////////////////////////////////
 ////  Boosting if needed   ///// 
////////////////////////////////
  
  
  double vx = P.get_x() / E;
  double vy = P.get_y() / E;
  double vz = P.get_z() / E;
  
  
  if ( std::abs(E*E - s) > 0.000000000001 )
  {
        
    CMatrix B = Boost(vx,vy,vz);
    p1 = B * p1;
    p2 = B * p2;
    
  }

  J = J*p;
  J = J/(4*PI*sqrt(s));

}

void Recursive_PSP(Vector P, int NPar, Vector* PList, double* MList, double* r, double& J){
    
   
    CMatrix Boosts[NPar];
    for(int i=0;i<NPar;i++){Boosts[i] = 1.;}
    
    
    int count = NPar;
    Vector P_Left = P;
    while (count > 2){
    
//     cout << "Entering recursive generation cycle, "<< count << " particles remaining" <<endl;
    double s = P_Left*P_Left;
    double E = P_Left.get_t();
    
//     cout << "Momemtum left = " << P_Left << endl;
    
    double s_min = 0;
    for(int i=0;i<(count-1);i++){ s_min += MList[i]; }
    s_min = s_min*s_min;
    double s_max = (sqrt(s) - MList[count-1])*(sqrt(s) - MList[count-1]);
    
    double s_now = (s_max-s_min)*r[3*(NPar-count)] + s_min;
    
    double E_now = (s - s_now + MList[count-1]*MList[count-1]) / (2*sqrt(s)) ;
    double pp_now = sqrt( E_now*E_now - MList[count-1]*MList[count-1] );
    double cos_th_now = 2 * r[3*(NPar-count)+1] - 1  ;
    double sin_th_now = sqrt( 1.0 - cos_th_now*cos_th_now );
    double ph_now = 2*PI*r[3*(NPar-count)+2] ;
    
    J = J*pp_now;
    J = J*(s_max-s_min);
    J = J/(8*PI*PI*sqrt(s));
    
    PList[count-1].set_t( E_now );
    PList[count-1].set_x( pp_now*cos(ph_now)*sin_th_now );
    PList[count-1].set_y( pp_now*sin(ph_now)*sin_th_now );
    PList[count-1].set_z( pp_now*cos_th_now );
    
//     cout << "p" << count << " generated!" << endl;
//     cout << "p" << 2*NPar-count << " =" << PList[2*NPar-count-1] << endl;
    
    CMatrix Boost_now;
    Boost_now = 1.;
    
    if ( std::abs(E*E - s) > 0.000000000001 ){
    
    double vx = P_Left.get_x() / E;
    double vy = P_Left.get_y() / E;
    double vz = P_Left.get_z() / E;  
//     cout << "v =("<<vx<<","<<vy<<","<<vz<<")"<<endl;
    Boost_now = Boost(vx,vy,vz);
    
    }
    
//     cout << "Cumulative Boosts" <<endl;
    for(int i=0;i<=(count-1);i++){
//         cout << "Boost[" << i << "] =" << endl << Boosts[i] << endl;
//         cout << "Boost_now =" << endl << Boost_now << endl;
        Boosts[i] = Boosts[i]*Boost_now;
//         cout << "Boost[" << i << "] =" << endl << Boosts[i] << endl;
    }
    
    P_Left = sqrt(s)*P0 - PList[count-1]; 
    
    
    count = count - 1;
    }
    
    
    double ran[2] = {r[3*NPar-6],r[3*NPar-5]};
    Phase_Space_Point_Generator_2_Particle(P_Left,PList[0],MList[0],PList[1],MList[1],ran,J);
    
    for(int i=0;i<NPar;i++){PList[i] = Boosts[i]*PList[i];}
    
    
    
    
}


    /////////////////////////////////////////////////////////////
   ///                                                       ///
  ///       Dedicated Dipole Phase Space Generator          ///
 ///                                                       ///
/////////////////////////////////////////////////////////////


double sig_of_min(double s, double Ma, double Mb, double Mi, double Mia){
    
    double aux[5];
    
    aux[0] = (s + Ma*Ma - Mb*Mb) / (2*sqrt(s));
    aux[1] = (1./2)*sqrt(lambda(s,Ma*Ma,Mb*Mb)/s);
    
    aux[2] = (aux[0]*aux[0]*(s + Mi*Mi - Mia*Mia)*(s + Mi*Mi - Mia*Mia)) - (aux[1]*aux[1]*(s - Mi*Mi + Mia*Mia)*(s - Mi*Mi + Mia*Mia));
    
    aux[3] = s*(aux[0]*aux[0]*(s + Mi*Mi - Mia*Mia)*(s + Mi*Mi - Mia*Mia) - aux[1]*aux[1]*lambda(s,Mi*Mi,Mia*Mia));
    aux[3] = 2*aux[0]*Mia*(s + Mi*Mi - Mia*Mia)*sqrt(aux[3]);
    
    aux[4] = Ma*Ma*(s+Mia*Mia)*(s + Mi*Mi - Mia*Mia)*(s + Mi*Mi - Mia*Mia) + 4*s*s*aux[1]*aux[1]*Mi*Mi;
    
    aux[4] = aux[4] - aux[3] ;
    
    aux[4] = aux[4] / aux[2] ; 
    
    return aux[4];    
    
    
    
}


double xia_min(double s, double Ma, double Mb, double Mi, double Mia){
    double aux[5] ; 
    
    aux[0] = (s + Ma*Ma - Mb*Mb)/(2*sqrt(s));
    aux[1] = (1./2)*sqrt(lambda(s,Ma*Ma,Mb*Mb)/s);
    
    aux[2] = sig_of_min(s,Ma,Mb,Mi,Mia);
    
    aux[3] = aux[0]*(s + aux[2] - Mia*Mia) - aux[1]*sqrt(lambda(s,aux[2],Mia*Mia));
    
    aux[4] = sqrt(s)*(aux[2] - Mi*Mi);
    
    aux[4] = aux[4] / aux[3];
    
    aux[4] = 1.0 - aux[4];
    
    return aux[4];
    
    
    
    
    
    
    
    
    
    
    
}


double sig_bound(double s, double Ma, double Mb, double Mi, double Mia, double x, int pm){
    double aux[6];
    double pm_c = double(pm);
    
    aux[0] = (1.-x)*(s + Ma*Ma - Mb*Mb)/(2*sqrt(s));
    aux[1] = (1./2)*(1.-x)*sqrt(lambda(s,Ma*Ma,Mb*Mb)/s);
    
    aux[2] = aux[0] - sqrt(s);
    aux[3] = aux[0]*(s - Mia*Mia) + sqrt(s)*Mi*Mi;
    
    aux[4] = aux[2]*aux[2] - aux[1]*aux[1];
    aux[5] = 2*(aux[2]*aux[3] + aux[1]*aux[1]*(s + Mia*Mia));
    aux[0] = aux[3]*aux[3] - aux[1]*aux[1]*(s-Mia*Mia)*(s-Mia*Mia);
    
    aux[5] = (-aux[5] + pm_c*sqrt(aux[5]*aux[5] - 4*aux[4]*aux[0]))/(2*aux[4]);
    
    return aux[5];
    
    
    
    
}








#endif