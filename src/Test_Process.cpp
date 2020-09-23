
#include "nlox_olp.h"
#include "Phase_Space_Tools.h"
#include "processconst.h"

int main(){
    
    NLOX_OLP_Start(NULL,NULL);

    ProcessConst pc;
    
    int next = 4;
    FourVector P(1000,0,0,0);
    FourVector PIn[2];
    double MIn[2] = {pc.mb.real(),0};
    double rIn[2] = {0.0,0.0};
    FourVector POut[3];
    double MOut[3] = {pc.mZ.real(),pc.mb.real()};
    double rOut[3*next-10] = {0.932123743289,0.774327483269};
    double JIn=1,JOut=1;
    Recursive_PSP(P,2,PIn,MIn,rIn,JIn);
    Recursive_PSP(P,next-2,POut,MOut,rOut,JOut);
    
    FourVector p1 = PIn[0];
    FourVector p2 = PIn[1];
    FourVector p3 = POut[0];
    FourVector p4 = POut[1];
    
    std::vector<FourVector> p;
    p.push_back(p1);
    p.push_back(p2);
    p.push_back(p3);
    p.push_back(p4);

    double rval[2];
    double acc;

    char sub[] = "bg_Zb";
    char typ[] = "tree_tree";
    char cp[] = "as1ae1";
    double mu = 500.;

    NLOX_OLP_EvalSubProcess(sub,typ,cp,p,&next,&mu,rval,&acc);

    return 0;
}