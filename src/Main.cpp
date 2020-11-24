#include "Real.h"
#include "Virtual.h"
#include "Integrators.h"
#include "PDF_Sets.h"

#define GeV2toPB 3.89379338E+8

class XSection{

    public:

        Process * Proc;
        RealIntegrands * Reals;
        VirtualIntegrands * Virtuals;
        LHAPDF_Set * PDF;

        XSection(){
            PDF = new LHAPDF_Set("CT10nlo");
            Proc = new Process();
            Reals  = new RealIntegrands(Proc);
            Virtuals = new VirtualIntegrands(Proc);
        }

        ~XSection(){
            delete PDF;
            delete Virtuals;
            delete Reals;
            delete Proc;
        }

};

XSection * X1;
void XSectionStart(){
    X1 = new XSection();
}

double Integarnd(double* x,size_t dim, void* param){

    double sqrts = 13000;

    if(x[0]<1E-3||x[1]<1E-3)return 0;

    double xa,xb;
    // PDF Parametrization 
    // Here x[0] always corresponds to 

    double m1 = 4.2;
    double m2 = 0.0;
    double shatmin = (4.2+91.1876)*(4.2+91.1876);
    double sqrtshat = sqrt(m1*m1 + m2*m2 + m1*m1*m2*m2/(sqrts*sqrts*x[0]*x[1])+sqrts*sqrts*x[0]*x[1]);
    double Mu = X1->Proc->pc.mZ.real();
    double pdffac = (X1->PDF->Evaluate(0,x[0],Mu)*X1->PDF->Evaluate(5,x[1],Mu)+X1->PDF->Evaluate(5,x[0],Mu)*X1->PDF->Evaluate(0,x[1],Mu));
    
    // Matrix Element 

    double r[2] = {x[2],x[3]};
    double rval=1;
    if(sqrtshat<shatmin) return 0;
    X1->Virtuals->setECM(sqrtshat);
    X1->Virtuals->Born("bg_Zb","as1ae1",r,&rval);

    if(isnan(rval)){
        std::cout<<"NLOX has returned nan, backtracking..."<<std::endl;
        std::cout<<"r = {"<<r[0]<<","<<r[1]<<"}"<<std::endl;
        std::cout<<"sqrtshat = "<<sqrtshat<<std::endl;
    }

    // Flux Factor

    double flux = 1.0/(2*lambda(sqrtshat*sqrtshat,m1*m1,m2*m2));

    // std::cout 

    return pdffac*flux*rval*GeV2toPB;
}

int main(int argc, char* argv[]){

    XSectionStart();

    GSL_Integrator GIntegrator(Integarnd,4);
    montecarlo_specs mc;
    GIntegrator.Integrate(&mc,"Vegas");

    delete X1;
    
    return 0;
}