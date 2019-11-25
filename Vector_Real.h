#ifndef _Vector_Real_H
#define _Vector_Real_H


#include <iostream>
#include <cmath>

using namespace std;



///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////                                         ///////////////////////////
////////////////////////////              Lorentz Vectors            ///////////////////////////
///////////////////////////                                         ///////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////



class Vector{
    /// Components ///
    double t,x,y,z;
    /// Printing Function ///
    friend std::ostream &operator<<(std::ostream &, const Vector&);
    public:
        
        /// Constructors ///
        Vector();
        Vector(const long double ,const long double ,const long double ,const long double );
        Vector(const Vector&); /// Copy Constructor ///
        
        /// Assignment Operator /// 
        Vector& operator=(const Vector&);
        
        /// Deconstructor ///
        ~Vector();

        
        /// Member Functions ///
        
	// Cartesian //
	double get_t() const ;
        double get_x() const ;
        double get_y() const ;
        double get_z() const ;
        
	// Polar //
	double get_ph() const;
	double get_th() const;
	
	
        void set_t(const double);
        void set_x(const double);
        void set_y(const double);
        void set_z(const double);
        
                
        
        Vector operator + (const Vector&);
        Vector operator - (const Vector&);
        Vector operator!();
	Vector operator-();
};

Vector P0(1,0,0,0);
Vector P1(0,1,0,0);
Vector P2(0,0,1,0);
Vector P3(0,0,0,1);




/// Printing Function ///

std::ostream &operator<<(std::ostream &output, const Vector& c)
{
  double t = c.get_t();
  double x = c.get_x();
  double y = c.get_y();
  double z = c.get_z();
  output << "(" << t << ", " << x << " , " << y << " , " << z << ")";
  return output;  
}

Vector::Vector() 
{
    double aux=0;
    t = aux;
    x = aux;
    y = aux;
    z = aux;
}

Vector::Vector(const long double  tm, const long double  xm, const long double  ym, const long double  zm ) 
{
    double aux1 = tm;
    double aux2 = xm;
    double aux3 = ym;
    double aux4 = zm;
        
    set_t(aux1);
    set_x(aux2);
    set_y(aux3);
    set_z(aux4);
  
}


Vector::Vector(const Vector& c) 
{
  t = c.get_t();
  x = c.get_x();
  y = c.get_y();
  z = c.get_z();
}

Vector::~Vector() 
{
    double aux = 0;
    
    set_t(aux);
    set_x(aux);
    set_y(aux);
    set_z(aux);
}

Vector& Vector::operator=(const Vector& c) 
{
  t = c.get_t();
  x = c.get_x();
  y = c.get_y();
  z = c.get_z();
  return *this;
}

/// Member functions ///

double Vector::get_t() const {return t;}
double Vector::get_x() const {return x;}
double Vector::get_y() const {return y;}
double Vector::get_z() const {return z;}

double Vector::get_ph() const
{
 
  long double  re = x;
  long double  im = y;
  long double  phi = atan(im/re);
  long double  pi = M_PI;
  if (re == 0 && im == 0) {return 0;}
  if ( re >= 0 && im >= 0 ) {return phi;}
  if ( re >= 0 && im <= 0 ) {return 2*pi + phi;}
  if ( re <= 0 ) {return pi + phi;}
  
}

double Vector::get_th() const 
{
  float pres = 0.0000000000000001;
  float ax,ay,az;
  ax = x;
  ay = y;
  az = z;
  if ( ( std::abs(ax) < pres ) && ( std::abs(ay) < pres ) && ( std::abs(az) < pres ) ) return 0;
  else return acos( az / sqrt( (ax*ax) + (ay*ay) + (az*az) ) );
  
  
}

void Vector::set_t(const double tm) {t = tm;}
void Vector::set_x(const double xm) {x = xm;}
void Vector::set_y(const double ym) {y = ym;}
void Vector::set_z(const double zm) {z = zm;}


/////////////////////////////////////////////////////////////////////////////////////
//////////////////       Overloading of Operators         ////////////////////////// 
///////////////////////////////////////////////////////////////////////////////////



/////////////////////// Sum //////////////////////////////////////////

Vector Vector::operator + (const Vector& param)
{
    Vector sum;
    sum.t = t + param.t;
    sum.x = x + param.x;
    sum.y = y + param.y;
    sum.z = z + param.z;
    return sum;
}

//////////////////////// Subtraction /////////////////////////////////////

Vector Vector::operator - (const Vector& param)
{
    Vector sub;
    sub.t = t - param.t;
    sub.x = x - param.x;
    sub.y = y - param.y;
    sub.z = z - param.z;
    return sub;
}

//////////////////////////// Dot Product ///////////////////////////////////

double operator * (Vector& lhs, Vector& rhs )
{
    double prod = (lhs.get_t())*rhs.get_t();
    prod = prod - (lhs.get_x())*rhs.get_x();
    prod = prod - (lhs.get_y())*rhs.get_y();
    prod = prod - (lhs.get_z())*rhs.get_z();
    return prod;
}

/////////////////////////// Scalar Product  ////////////////////////////////

Vector operator * (long double  lhs, Vector& rhs)
{
    Vector prod(lhs*rhs.get_t(),lhs*rhs.get_x(),lhs*rhs.get_y(),lhs*rhs.get_z());
    return prod;  
    
}

Vector operator * (Vector& lhs, long double  rhs)
{
    Vector prod(rhs*lhs.get_t(),rhs*lhs.get_x(),rhs*lhs.get_y(),rhs*lhs.get_z());
    return prod;  
    
}

Vector operator * (double lhs, Vector& rhs)
{
    Vector prod(lhs*rhs.get_t(),lhs*rhs.get_x(),lhs*rhs.get_y(),lhs*rhs.get_z());
    return prod;  
    
}

Vector operator * (Vector& lhs, double rhs)
{
    Vector prod(rhs*lhs.get_t(),rhs*lhs.get_x(),rhs*lhs.get_y(),rhs*lhs.get_z());
    return prod;  
    
}


/////////////////////////// Comparison //////////////////////////////////////

bool operator == (Vector& lhs, Vector& rhs)
{
    return lhs.get_t()==rhs.get_t() && lhs.get_x()==rhs.get_x() && lhs.get_y()==rhs.get_y() && lhs.get_z()==rhs.get_z() ;
}

bool operator != (Vector& lhs, Vector& rhs)
{
    return 1-(rhs==lhs);
}

//////////////////////// Negative ///////////////////////////////////////////

Vector Vector::operator!()
{
    Vector cong(!t,!x,!y,!z);
    return cong;
    
    
}


Vector Vector::operator-()
{
  Vector neg(-t,-x,-y,-z);
  return neg;
  
}


//////////////////////////////////////////////////////////////////////////////
///////////////////////  Some Extra Functions  //////////////////////////////
////////////////////////////////////////////////////////////////////////////


int look_for(int b,int *a)
{
    int pos;
    int cont = 0;
    int i;
    for(i=0; i<4; i++)
        {
        if (b == a[i])
            {
            cont ++;
            pos = i;                
            }           
        }
    if (cont == 1) return pos;
    else return -1;
    
    
}

int sig(int *a)
{
    int signo = 1;
    int i;
    int vec[4];
    for (i=0; i<4; i++) vec[i] = a[i];
    for (i=0; i<4; i++)
        {
        if(vec[i] == i);
        else 
            {
            int pos = look_for(i,vec); 
            if( pos == -1 ) 
                {
                signo = 0;
                }
             else 
                {
                vec[pos] = vec[i];
                vec[i] = i;
                signo = (-1)*signo;
                }
             }
          }
    return signo;
    
}

int signo(int i,int j, int k, int l)
    {
    int a[4];
    a[0] = i;
    a[1] = j;
    a[2] = k;
    a[3] = l;
    return sig(a);
    }

double epsilon(Vector a, Vector b, Vector c, Vector d)
{
    double det = 0;
    double M[4][4];
    M[0][0] = a.get_t();
    M[0][1] = a.get_x();
    M[0][2] = a.get_y();
    M[0][3] = a.get_z();
    
    M[1][0] = b.get_t();
    M[1][1] = b.get_x();
    M[1][2] = b.get_y();
    M[1][3] = b.get_z();
    
    M[2][0] = c.get_t();
    M[2][1] = c.get_x();
    M[2][2] = c.get_y();
    M[2][3] = c.get_z();
    
    M[3][0] = d.get_t();
    M[3][1] = d.get_x();
    M[3][2] = d.get_y();
    M[3][3] = d.get_z();
    
    int i,j,k,l;
    for (i=0; i<4; i++)
        {
        for(j=0; j<4; j++)
            {
            for(k=0; k<4; k++)
                {
                for(l=0; l<4; l++)
                    {
                    det = det + (signo(i,j,k,l)*M[0][i]*M[1][j]*M[2][k]*M[3][l]);                       
                    }
                    
                }
                
                
            }
            
            
        }
    return det;
}




///////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////                                         ///////////////////////////
////////////////////////////                 Matrices                ///////////////////////////
///////////////////////////                                         ///////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////


class CMatrix
{
    
    
    /// Printing Function ///
    friend std::ostream &operator<<(std::ostream &, const CMatrix&);
    
    public:
      
       /// Components ///
       double M[4][4];
       
       
       /// Constructors ///
       CMatrix();
       CMatrix(const CMatrix&); /// Copy Constructor ///  

       /// Assignment Operator ///
       CMatrix& operator=(const CMatrix&);
       CMatrix& operator=(const double& );
       
       /// Deconstructor ///
       ~CMatrix();       

       /// Overloaded Operators ///
       CMatrix operator+(const CMatrix&) const;
       CMatrix operator-(const CMatrix&) const;
       CMatrix operator*(const CMatrix&);
       
       


};


///  Printing Function   ///

std::ostream &operator<<(std::ostream &output, const CMatrix& c)
{
    int i,j;
    for(j=0; j<4; j++)
    {
        output << "| ";
        for(i=0; i<4; i++)
        {
            output << c.M[j][i] << " ";
        
        }
        output << "|\n";
    }
    return output;
    
    
}


///  Constructors ///

CMatrix::CMatrix()
{
    int i,j;
    for(i=0; i<4; i++)
    {
        for(j=0; j<4; j++)
        {
            M[i][j] = 0.;
            
        }
              
    }  
    
}

CMatrix::~CMatrix()
{
  int i,j;
    for(i=0; i<4; i++)
    {
        for(j=0; j<4; j++)
        {
            M[i][j] = 0.;
            
        }
              
    }  
      
    
    
    
}

CMatrix::CMatrix(const CMatrix& c)
{
  int i,j;
    for(i=0; i<4; i++)
    {
        for(j=0; j<4; j++)
        {
            M[i][j] = c.M[i][j];
            
        }
              
    }  
      
    
    
    
}


CMatrix& CMatrix::operator=(const CMatrix& c) 
{
  int i,j;
    for(i=0; i<4; i++)
    {
        for(j=0; j<4; j++)
        {
            M[i][j] = c.M[i][j];
            
        }
              
    }   
  return *this;
}


CMatrix& CMatrix::operator=(const double & c) 
{
    double aux = c;
    int i,j;
    for(i=0; i<4; i++)
    {
        M[i][i] = aux;
        
    }   
    return *this;
}

/////////////////////////////////////////////////////////////////////////////////////
//////////////////       Overloading of Operators         ////////////////////////// 
///////////////////////////////////////////////////////////////////////////////////


CMatrix CMatrix::operator+(const CMatrix& c) const
{
    CMatrix sum;
    int i,j;
    for(i=0; i<4; i++)
    {
        for(j=0; j<4; j++)
        {
            double entry = M[i][j] + c.M[i][j];
            sum.M[i][j] = entry;
                       
        }
              
    }  
   return sum; 
    
    
    
}

CMatrix CMatrix::operator-(const CMatrix& c) const
{
   CMatrix sub;
    int i,j;
    for(i=0; i<4; i++)
    {
        for(j=0; j<4; j++)
        {
            double entry = M[i][j] - c.M[i][j];
            sub.M[i][j] = entry;
                       
        }
              
    }  
   return sub; 
     
    
    
    
    
    
    
}

CMatrix CMatrix::operator*(const CMatrix& c)
{
    CMatrix prod;
    int i,j,k,l;
    for(i=0; i<4; i++)
    {
        for(j=0; j<4; j++)
        {
            double entry=0.;
            for(l=0; l<4; l++)
            {
                double a = c.M[l][j];
                double b = M[i][l]*a;
                entry = entry + b; 
                
            }  
            prod.M[i][j] = entry;
                       
        }
              
    }  
   return prod;     
    
}


///  Useful operators ///

Vector operator * (CMatrix lhs, Vector rhs)
{
  Vector aux;
  aux.set_t( lhs.M[0][0]*rhs.get_t()+lhs.M[0][1]*rhs.get_x()+lhs.M[0][2]*rhs.get_y()+lhs.M[0][3]*rhs.get_z() );
  aux.set_x( lhs.M[1][0]*rhs.get_t()+lhs.M[1][1]*rhs.get_x()+lhs.M[1][2]*rhs.get_y()+lhs.M[1][3]*rhs.get_z() );
  aux.set_y( lhs.M[2][0]*rhs.get_t()+lhs.M[2][1]*rhs.get_x()+lhs.M[2][2]*rhs.get_y()+lhs.M[2][3]*rhs.get_z() );
  aux.set_z( lhs.M[3][0]*rhs.get_t()+lhs.M[3][1]*rhs.get_x()+lhs.M[3][2]*rhs.get_y()+lhs.M[3][3]*rhs.get_z() );
  return aux; 
  
}


CMatrix operator * (double c, CMatrix m)
{
    CMatrix out;
    int i,j;
    for(i=0; i<4; i++)
    {
        for(j=0; j<4; j++)
        {
            out.M[i][j] = c*m.M[i][j];        
            
        }
        
        
    }
    
   return out; 
    
}

CMatrix operator * (long double  c, CMatrix m)
{
    CMatrix out;
    int i,j;
    for(i=0; i<4; i++)
    {
        for(j=0; j<4; j++)
        {
            out.M[i][j] = c*m.M[i][j];        
            
        }
        
        
    }
    
   return out; 
    
}

CMatrix operator * (CMatrix m, double c)
{
    CMatrix out;
    int i,j;
    for(i=0; i<4; i++)
    {
        for(j=0; j<4; j++)
        {
            out.M[i][j] = c*m.M[i][j];        
            
        }
        
        
    }
    
   return out; 
    
}

CMatrix operator * (CMatrix m, long double  c)
{
    CMatrix out;
    int i,j;
    for(i=0; i<4; i++)
    {
        for(j=0; j<4; j++)
        {
            out.M[i][j] = c*m.M[i][j];        
            
        }
        
        
    }
    
   return out; 
    
}


CMatrix operator + (double c, CMatrix m)
{
    CMatrix out;
    int i,j;
    for(i=0; i<4; i++)
    {
        for(j=0; j<4; j++)
        {
            out.M[i][j] = m.M[i][j];
            if (i == j ) out.M[i][i] = out.M[i][i] + c;
            
        }
        
        
    }
    
   return out; 
    
}

CMatrix operator + (long double  c, CMatrix m)
{
    CMatrix out;
    int i,j;
    for(i=0; i<4; i++)
    {
        for(j=0; j<4; j++)
        {
            out.M[i][j] = m.M[i][j];
            if (i == j ) out.M[i][i] = out.M[i][i] + c;
            
        }
        
        
    }
    
   return out; 
    
}

CMatrix operator + (CMatrix m, double c)
{
    CMatrix out;
    int i,j;
    for(i=0; i<4; i++)
    {
        for(j=0; j<4; j++)
        {
            out.M[i][j] = m.M[i][j];
            if (i == j ) out.M[i][i] = out.M[i][i] + c;
            
        }
        
        
    }
    
   return out; 
    
}

CMatrix operator + (CMatrix m, long double  c)
{
    CMatrix out;
    int i,j;
    for(i=0; i<4; i++)
    {
        for(j=0; j<4; j++)
        {
            out.M[i][j] = m.M[i][j];
            if (i == j ) out.M[i][i] = out.M[i][i] + c;
            
        }
        
        
    }
    
   return out; 
    
}

CMatrix operator - (double c, CMatrix m)
{
    CMatrix out;
    int i,j;
    for(i=0; i<4; i++)
    {
        for(j=0; j<4; j++)
        {
            out.M[i][j] = -m.M[i][j];
            if (i == j ) out.M[i][i] = out.M[i][i] + c;
            
        }
        
        
    }
    
   return out; 
    
}

CMatrix operator - (long double  c, CMatrix m)
{
    CMatrix out;
    int i,j;
    for(i=0; i<4; i++)
    {
        for(j=0; j<4; j++)
        {
            out.M[i][j] = -m.M[i][j];
            if (i == j ) out.M[i][i] = out.M[i][i] + c;
            
        }
        
        
    }
    
   return out; 
    
}

CMatrix operator - (CMatrix m, double c)
{
    CMatrix out;
    int i,j;
    for(i=0; i<4; i++)
    {
        for(j=0; j<4; j++)
        {
            out.M[i][j] = m.M[i][j];
            if (i == j ) out.M[i][i] = out.M[i][i] - c;
            
        }
        
        
    }
    
   return out; 
    
}

CMatrix operator - (CMatrix m, long double  c)
{
    CMatrix out;
    int i,j;
    for(i=0; i<4; i++)
    {
        for(j=0; j<4; j++)
        {
            out.M[i][j] = m.M[i][j];
            if (i == j ) out.M[i][i] = out.M[i][i] - c;
            
        }
        
        
    }
    
   return out; 
    
}

CMatrix operator / (CMatrix m, double c)
{
    CMatrix out;
    int i,j;
    for(i=0; i<4; i++)
    {
        for(j=0; j<4; j++)
        {
            out.M[i][j] = m.M[i][j]/c;
            
        }
        
        
    }
    
   return out; 
    
}

CMatrix operator / (CMatrix m, long double  c)
{
    CMatrix out;
    int i,j;
    for(i=0; i<4; i++)
    {
        for(j=0; j<4; j++)
        {
            out.M[i][j] = m.M[i][j]/c;
            
        }
        
        
    }
    
   return out; 
    
}




//////////////////////////////////////////////////////////////////////////////
///////////////////////  Some Extra Functions  //////////////////////////////
////////////////////////////////////////////////////////////////////////////



double Tr(CMatrix c)
{
    double trace;
    int i,j;
    for(i=0; i<4; i++)
    {
    trace = trace + c.M[i][i];
              
    }  
    return trace;
    
    
}

double Det(CMatrix c)
{
    double det;
    int i,j,k,l;
    for (i=0; i<4; i++)
        {
        for(j=0; j<4; j++)
            {
            for(k=0; k<4; k++)
                {
                for(l=0; l<4; l++)
                    {
                    double sig = signo(i,j,k,l);
                    det = det + (sig*c.M[0][i]*c.M[1][j]*c.M[2][k]*c.M[3][l]);                       
                    }
                    
                }
                
                
            }
            
            
        }
    return det;
    
}









#endif
