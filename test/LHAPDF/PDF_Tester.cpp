#include "../../src/PDF_Sets.h"

int main(){

    LHAPDF_Set MYPDF("cteq6");

    int PID = 0; // Gluon
    double MuFact = 100; // ~mZ

    for(double x=1E-3;x<1;x+=1E-3){
        std::cout<<MYPDF.Evaluate(PID,x,MuFact)<<std::endl;
    }


  return 0;
}