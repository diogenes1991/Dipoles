#ifndef F_VECTOR_H_
#define F_VECTOR_H_

#include <iostream>
#include <sstream>
#include <string>
#include <iomanip>
#include <cmath>

/**********************************************************************
 *                                                                    *
 * FVectorT                                                         *
 *                                                                    *
 * Simple implementation of Minkowski vectors, including basic        *
 * algebra and publicly accessible members.                           *
 *                                                                    *
 **********************************************************************/

template <typename fT> class FVectorT {
  public:
    FVectorT() : p0(0), p1(0), p2(0), p3(0) {}

    FVectorT(const FVectorT & fv)
    : p0(fv.p0), p1(fv.p1), p2(fv.p2), p3(fv.p3) {}

    template <typename T2> FVectorT<fT>(const FVectorT<T2> & fv) {
      p0 = fT(fv.p0);
      p1 = fT(fv.p1);
      p2 = fT(fv.p2);
      p3 = fT(fv.p3);
    }

    FVectorT(fT x0, fT x1, fT x2, fT x3)
    : p0(x0), p1(x1), p2(x2), p3(x3) {}

    void boost(int dir, fT gamma, fT gammabeta);
    void rotate_xy(fT vz, fT sint, fT cost, fT sinp, fT cosp);

    // operator overloading for standard algebra
    FVectorT<fT>  operator+(const FVectorT<fT> & fv);
    FVectorT<fT>& operator+=(const FVectorT<fT> & fv);
    FVectorT<fT>  operator-(const FVectorT<fT> & fv);
    FVectorT<fT>& operator-=(const FVectorT<fT> & fv);

    // Negation operator
    FVectorT<fT> operator-() const;

    FVectorT<fT> negate(); // Redundant
    FVectorT<fT> negateSpatial();

    fT comp(int i);
    fT abs();

    // Computes the scalar product between two FVectorTs.
    // ( signature g_mu,nu = (+1,-1,-1,-1) )
    // TODO: put outside class?
    static fT scalar_product(const FVectorT<fT> & v1,
                             const FVectorT<fT> & v2);

    // Computes the contraction of F FVectorTs with the
    // totally antisymmetric Levi-Civita tensor e_(mu,nu,ro.sig)
    static fT antisym_product(const FVectorT<fT> & v1,
                              const FVectorT<fT> & v2,
		              const FVectorT<fT> & v3,
		              const FVectorT<fT> & v4);

    void print(std::ostream & strm = std::cout, std::string prefix = "");
    void print(int prec, int width, std::ostream & strm = std::cout,
               std::string prefix = "");

    std::string recolaString(std::string suffix = "_dp",
                             bool split = false, std::string padding = "");

    // public, we know what we do
    fT p0, p1, p2, p3;
};

typedef FVectorT<double> FVector;

// Template functions

template <typename fT> inline FVectorT<fT> FVectorT<fT>::operator+(
    const FVectorT<fT> & fv) {
  return FVectorT<fT>(p0 +fv.p0, p1 + fv.p1, p2 + fv.p2, p3 + fv.p3);
}

template <typename fT> inline FVectorT<fT>& FVectorT<fT>::operator+=(
    const FVectorT<fT> & fv) {
  p0 += fv.p0;
  p1 += fv.p1;
  p2 += fv.p2;
  p3 += fv.p3;
  return *this;
}

template <typename fT> inline FVectorT<fT> FVectorT<fT>::operator-(
    const FVectorT<fT> & fv) {
  return FVectorT<fT>(p0 - fv.p0, p1 - fv.p1, p2 - fv.p2, p3 - fv.p3);
}

template <typename fT> inline FVectorT<fT>& FVectorT<fT>::operator-=(
    const FVectorT<fT> & fv) {
  p0 -= fv.p0;
  p1 -= fv.p1;
  p2 -= fv.p2;
  p3 -= fv.p3;
  return *this;
}

template <typename fT> inline FVectorT<fT>
FVectorT<fT>::operator-() const {
  return FVectorT<fT>(-p0, -p1, -p2, -p3);
}

template <typename fT> inline FVectorT<fT> operator*(
    const FVectorT<fT> & vec, const fT & val) {
  return FVectorT<fT>(val*vec.p0, val*vec.p1, val*vec.p2, val*vec.p3);
}

template <typename fT> inline FVectorT<fT> FVectorT<fT>::negate() {
  return FVectorT<fT>(-p0, -p1, -p2, -p3);
}

template <typename fT> inline FVectorT<fT>
FVectorT<fT>::negateSpatial() {
  return FVectorT<fT>(p0, -p1, -p2, -p3);
}

template <typename fT> inline fT FVectorT<fT>::abs() {
  return p0*p0 - p1*p1 - p2*p2 - p3*p3;
}

template <typename fT> inline fT FVectorT<fT>::comp(int i) {
  fT retval = 0;
  switch (i) {
    case 0: retval = p0;
      break;
    case 1: retval = p1;
      break;
    case 2: retval = p2;
      break;
    case 3: retval = p3;
      break;
    default: std::cerr << "error: FVectorT::comp: index " << i
                       << " ne [0,3]" << std::endl;
      exit(1);
    }
  return retval;
}

template <typename fT> inline fT FVectorT<fT>::scalar_product(
    const FVectorT & v1, const FVectorT & v2) {
  return v1.p0*v2.p0 - v1.p1*v2.p1 - v1.p2*v2.p2 - v1.p3*v2.p3;
}

// non-inlined functions

template <typename fT> fT FVectorT<fT>::antisym_product(
    const FVectorT & v1, const FVectorT & v2,
    const FVectorT & v3, const FVectorT & v4) {
  return + v1.p0*v2.p1 * (v3.p2*v4.p3 - v3.p3*v4.p2)
         + v1.p0*v2.p2 * (v3.p3*v4.p1 - v3.p1*v4.p3)
         + v1.p0*v2.p3 * (v3.p1*v4.p2 - v3.p2*v4.p1)
         + v1.p1*v2.p0 * (v3.p3*v4.p2 - v3.p2*v4.p3)
         + v1.p1*v2.p2 * (v3.p0*v4.p3 - v3.p3*v4.p0)
         + v1.p1*v2.p3 * (v3.p2*v4.p0 - v3.p0*v4.p2)
         + v1.p2*v2.p0 * (v3.p1*v4.p3 - v3.p3*v4.p1)
         + v1.p2*v2.p1 * (v3.p3*v4.p0 - v3.p0*v4.p3)
         + v1.p2*v2.p3 * (v3.p0*v4.p1 - v3.p1*v4.p0)
         + v1.p3*v2.p0 * (v3.p2*v4.p1 - v3.p1*v4.p2)
         + v1.p3*v2.p1 * (v3.p0*v4.p2 - v3.p2*v4.p0)
         + v1.p3*v2.p2 * (v3.p1*v4.p0 - v3.p0*v4.p1);
}

template <typename fT> void FVectorT<fT>::boost(int dir,
                                                   fT gamma, fT gammabeta) {
  fT pp0, pp1, pp2, pp3;
  switch(dir) {
    case 1: pp0 = gamma*p0 + gammabeta*p1;
            pp1 = gammabeta*p0 + gamma*p1;
            pp2 = p2;
            pp3 = p3;
            break;
    case 2: pp0 = gamma*p0 + gammabeta*p2;
            pp1 = p2;
            pp2 = gammabeta*p0 + gamma*p2;
            pp3 = p3;
            break;
    case 3: pp0 = gamma*p0 + gammabeta*p3;
            pp1 = p1;
            pp2 = p2;
            pp3 = gammabeta*p0 + gamma*p3;
            break;
  }
  p0 = pp0;
  p1 = pp1;
  p2 = pp2;
  p3 = pp3;
}

template <typename fT> void FVectorT<fT>::rotate_xy(fT vz,
                                                       fT sint, fT cost,
                                                       fT sinp, fT cosp) {
  fT pp1, pp2, pp3;
  pp1 = vz*(cosp*(cost*p1 + sint*p3) - sinp*p2);
  pp2 = vz*(sinp*(cost*p1 + sint*p3) + cosp*p2);
  pp3 = vz*(-sint*p1 + cost*p3);

  p1 = pp1;
  p2 = pp2;
  p3 = pp3;
}

template <typename fT> void FVectorT<fT>::print(std::ostream & strm,
                                                   std::string prefix) {
  strm << prefix << p0 << " " << p1 << " " << p2 << " " << p3 << std::endl;
}

template <typename fT> void FVectorT<fT>::print(int prec, int width,
                                                   std::ostream & strm,
                                                   std::string prefix) {
  strm << std::setprecision(prec)
       << prefix
       << std::setw(width) << p0 << "\t"
       << std::setw(width) << p1 << "\t"
       << std::setw(width) << p2 << "\t"
       << std::setw(width) << p3 << std::endl;
}

///
///  Extended operators
///

template <class T> 
inline T operator * (FVectorT<T>& lhs, FVectorT<T>& rhs ){
    T prod = (lhs.p0)*rhs.p0;
    prod = prod - (lhs.p1)*rhs.p1;
    prod = prod - (lhs.p2)*rhs.p2;
    prod = prod - (lhs.p3)*rhs.p3;
    return prod;
}

template <class T,class U> 
inline FVectorT<T> operator * (const U lhs, FVectorT<T>& rhs){
    FVectorT<T> prod(lhs*rhs.p0,lhs*rhs.p1,lhs*rhs.p2,lhs*rhs.p3);
    return prod;  
    
}

template <class T,class U>
inline FVectorT<T> operator * (FVectorT<T>& lhs, const U rhs){
    FVectorT<T> prod(rhs*lhs.p0,rhs*lhs.p1,rhs*lhs.p2,rhs*lhs.p3);
    return prod;
}

/// 
/// FMatricesT Operators
/// 

template <class T>
class FMatrixT{
        
    public:
       T M[4][4];
    
       /// Constructors ///
       FMatrixT<T>();
       FMatrixT<T>(const FMatrixT<T>&); /// Copy Constructor ///  

       /// Assignment Operator ///
       FMatrixT<T>& operator=(const FMatrixT<T>&);
       FMatrixT<T>& operator=(const T&);
       
       /// Deconstructor ///
       ~FMatrixT<T>();       

       /// Overloaded Operators ///
       FMatrixT<T> operator+(const FMatrixT<T>&) const;
       FMatrixT<T> operator-(const FMatrixT<T>&) const;
       FMatrixT<T> operator*(const FMatrixT<T>&);


};

typedef FMatrixT<double> FMatrix;

template <class T>
std::ostream &operator<<(std::ostream &output, const FMatrixT<T>& c){
    for(int j=0; j<4; j++){
        output << "| ";
        for(int i=0; i<4; i++){
            output << c.M[j][i] << " ";
        }
        output << "|\n";
    }
    return output;
}

template <class T>
FMatrixT<T>::FMatrixT(){
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            M[i][j] = 0.;
        }
    } 
}

template <class T>
FMatrixT<T>::~FMatrixT(){
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            M[i][j] = 0.;
        }
    } 
}

template <class T>
FMatrixT<T>::FMatrixT(const FMatrixT& c){
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            M[i][j] = c.M[i][j];
        }
    }  
}

template <class T>
FMatrixT<T>& FMatrixT<T>::operator=(const FMatrixT<T>& c){
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            M[i][j] = c.M[i][j];
        }
    }   
  return *this;
}

template <class T>
FMatrixT<T>& FMatrixT<T>::operator=(const T& c) {
    T aux = c;
    for(int i=0; i<4; i++){
        M[i][i] = aux;
    }   
    return *this;
}
 
template <class T>
FMatrixT<T> FMatrixT<T>::operator+(const FMatrixT<T>& c)const{
    FMatrixT<T> sum;
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            T entry = M[i][j] + c.M[i][j];
            sum.M[i][j] = entry;
        }
    }  
   return sum; 
}

template <class T>
FMatrixT<T> FMatrixT<T>::operator-(const FMatrixT<T>& c)const{
    FMatrixT<T> sub;
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            T entry = M[i][j] - c.M[i][j];
            sub.M[i][j] = entry;
        }
    }  
   return sub; 
}

template <class T> 
FMatrixT<T> FMatrixT<T>::operator*(const FMatrixT<T>& c){
    FMatrixT<T> prod;
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            T entry=0.;
            for(int l=0; l<4; l++){
                T a = c.M[l][j];
                T b = M[i][l]*a;
                entry = entry + b; 
            }  
            prod.M[i][j] = entry;
        }
    }  
   return prod;     
}

template <class T>
FVectorT<T> operator * (FMatrixT<T> lhs, FVectorT<T> rhs){
  FVectorT<T> aux(lhs.M[0][0]*rhs.p0+lhs.M[0][1]*rhs.p1+lhs.M[0][2]*rhs.p2+lhs.M[0][3]*rhs.p3,
                     lhs.M[1][0]*rhs.p0+lhs.M[1][1]*rhs.p1+lhs.M[1][2]*rhs.p2+lhs.M[1][3]*rhs.p3,
                     lhs.M[2][0]*rhs.p0+lhs.M[2][1]*rhs.p1+lhs.M[2][2]*rhs.p2+lhs.M[2][3]*rhs.p3,
                     lhs.M[3][0]*rhs.p0+lhs.M[3][1]*rhs.p1+lhs.M[3][2]*rhs.p2+lhs.M[3][3]*rhs.p3);
  return aux; 
}

template <class T,class U>
FMatrixT<T> operator * (U c, FMatrixT<T> m){
    FMatrixT<T> out;
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            out.M[i][j] = c*m.M[i][j];        
        }
    }
   return out;
}

template <class T, class U>
FMatrixT<T> operator * (FMatrixT<T> m, U c){
    FMatrixT<T> out;
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            out.M[i][j] = c*m.M[i][j];        
        }
    }
   return out; 
}

template <class T, class U>
FMatrixT<T> operator + (U c, FMatrixT<T> m){
    FMatrixT<T> out;
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            out.M[i][j] = m.M[i][j];
            if (i == j ) out.M[i][i] = out.M[i][i] + c;
        }
    }
   return out; 
}

template <class T, class U>
FMatrixT<T> operator + (FMatrixT<T> m, U c){
    FMatrixT<T> out;
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            out.M[i][j] = m.M[i][j];
            if (i == j ) out.M[i][i] = out.M[i][i] + c;
        }
    }
   return out; 
}

template <class T, class U>
FMatrixT<T> operator - (U c, FMatrixT<T> m){
    FMatrixT<T> out;
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            out.M[i][j] = -m.M[i][j];
            if (i == j ) out.M[i][i] = out.M[i][i] + c;
        }
    }
   return out; 
}

template <class T, class U>
FMatrixT<T> operator - (FMatrixT<T> m, U c){
    FMatrixT<T> out;
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            out.M[i][j] = m.M[i][j];
            if (i == j ) out.M[i][i] = out.M[i][i] - c;
        }
    }
   return out; 
}

template <class T, class U>
FMatrixT<T> operator / (FMatrixT<T> m, U c){
    FMatrixT<T> out;
    for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
            out.M[i][j] = m.M[i][j]/c;
        }
    }
   return out; 
}

template<class T>
FMatrixT<T> TensorProduct(FVectorT<T> a, FVectorT<T> b){
    FMatrixT<T> Mat;
    T flata[4] = {a.p0,a.p1,a.p2,a.p3};
    T flatb[4] = {b.p0,b.p1,b.p2,b.p3};
    for(int i=0;i<4;i++){
        for(int j=0;j<4;j++){
            Mat.M[i][j] = flata[i]*flatb[j];
        }
    }
    return Mat;
}

template <class T>
T Trace(FMatrixT<T> M){
    return M.M[0][0] + M.M[1][1] + M.M[2][2] + M.M[3][3];
}

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


#endif // FOUR_VECTOR_H_
