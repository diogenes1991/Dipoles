#include <iostream>
#include <ctime>
#include <cmath>
#include "/home/diogenes1991/NLOX/tests/AA_bbbar_NEW/AA_bbbar_EW/code/olp.h"

#define nPar 4

double BORN(int NPar, double* PLIST){
    std::vector<FourVector> P;
    for(int i=0;i<NPar-1;i++){
        FourVector aux(PLIST[4*i],PLIST[4*i+1],PLIST[4*i+2],PLIST[4*i+3]);
        P.push_back(aux);
    }
    proc->setPSP2(P);
    return proc->born(1.0);
}

double VIRT(int NPar, double* PLIST){
    std::vector<FourVector> P;
    for(int i=0;i<NPar-1;i++){
        FourVector aux(PLIST[4*i],PLIST[4*i+1],PLIST[4*i+2],PLIST[4*i+3]);
        P.push_back(aux);
    }
    proc->setPSP2(P);
    double virtc2,virtc1,virtc0;
    double CTc2,CTc1,CTc0;
    proc->virt(proc->pc.mu, &virtc2, &virtc1, &virtc0, -1.0);
    proc->CT(proc->pc.mu, &CTc2, &CTc1, &CTc0, 1.0);
    return virtc0+CTc0;
}

using namespace std;

int main(int argc, char *argv[]){
    
//     double MOMENTA[4*nPar];
    
    if(argc<12) {
        cout << "Error: Not enough arguments" <<endl;
        return 0;
    }
    else{
        for(int i=0;i<4*(nPar-1);i++){
//          MOMENTA[i] = atoi(argv[i+1]);   
        }
    }
    
    
    double MOMENTA[4*4] = {250, 0 , 0 , 250,250, 0 , 0 , -250,250, 202.6971329859771 , 52.31723006465615 , 136.5984616224757,250, -202.6971329859771 , -52.31723006465615 , -136.5984616224757};
    
//     FourVector p1(250, 0 , 0 , 250);
//     FourVector p2(250, 0 , 0 , -250);
//     FourVector p3(250, 202.6971329859771 , 52.31723006465615 , 136.5984616224757);
//     FourVector p4(250, -202.6971329859771 , -52.31723006465615 , -136.5984616224757);
  
    cout << BORN(nPar,MOMENTA) << endl;
//     cout << "Virtual = " << VIRT(nPar,MOMENTA) << endl;


  
  
}
