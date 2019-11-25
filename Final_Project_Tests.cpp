#include <iostream>
#include <fstream>
#include <cmath>
#include "Vector.h"
#include "Standard_Model.h"
#include "Phase_Spaces.h"


using namespace std;


int DEBUG = 1;


double Charge( int PID ){
    double aux = 0;
    double Qu = 2/3.;
    if ( abs(PID) <= 6 && int(std::abs(PID))%2==1 ) aux = Qu-1.;
    if ( std::abs(PID) <= 6 && int(std::abs(PID))%2==0 ) aux = Qu;
    if ( std::abs(PID) == 11 || std::abs(PID) == 13 || std::abs(PID) == 15 ) aux = -1;
    if ( std::abs(PID) == 24 ) aux = 1;
    if ( PID < 0 ) aux *= -1.0; 
    return aux;
    
}


bool IsColored(int PID){
    if(abs(PID)<=6||PID==21) return 1;
    else return 0;
}


int signo( int j ){
    if ( j >= 0 ) return 1;
    else return -1;
    
}

string PID_to_Name(int PID){
    string out;
    int aux = abs(PID);
    /// Quarks ///
    if(aux<=6){
    if (aux==1) out += "d";
    if (aux==2) out += "u";
    if (aux==3) out += "s";
    if (aux==4) out += "c";
    if (aux==5) out += "b";
    if (aux==6) out += "t";
    if (PID<0) out += "bar";
    }
    /// Leptons ///
    else if(aux>10 && aux<17){
    if(aux%2){
        if(aux==11) out += "e";
        if(aux==13) out += "mu";
        if(aux==15) out += "tau";
        if(PID>0) out += "m";
        else out += "p";
    }
    else{
        out += "nu";
        if(aux==12) out += "e";
        if(aux==14) out += "mu";
        if(aux==16) out += "tau";
        if(PID<0) out += "bar";
    }
    }
    /// Gauge Bosons and Higgs ///
    else if( aux>20 &&  aux<26 ){
        if(aux==21) out+= "g";
        if(aux==22) out+= "A";
        if(aux==23) out+= "Z";
        if(PID==24) out+= "Wp";
        if(PID==-24) out+= "Wp";
        if(aux==25) out+= "H";
        
    }
    return out;
    
}


void Create_Subtracted(int NEXT, int* EXTID, int* RETCONFIG){
    
    if(DEBUG){
        for(int i=0;i<=NEXT-1;i++){
    if(i==0) cout << "QED dipole routine triggered for the process " << endl;
    if(i==2) cout << " ---> ";
    cout << PID_to_Name(EXTID[i]) << "(" << i+1 << ")";
    if(i!=NEXT-1&&i!=1) cout << " + ";
    else if(i!=1) cout << endl;
    }
    }

    double aux[2];
    
    bool USEALTCH = 0;
    
    for(int i=0;i<=NEXT-1;i++){
        for(int j=0;j<=NEXT-1;j++){
            double ch;
            if(USEALTCH) ch = Charge(EXTID[i])*Charge(EXTID[j])*signo(EXTID[i])*signo(EXTID[j]);
            else ch = Charge(abs(EXTID[i]))*Charge(abs(EXTID[j]))*signo(EXTID[i])*signo(EXTID[j]);
            if(ch&&j!=i){
//             double pe[5],ps[5];
//             for(int k=0;k<5;k++){pe[k]=pp[5*i+k];ps[k]=pp[5*j+k];}
            int conf = 0;
            if(i<2&&j<2) conf = 0;
            if(i<2&&j>=2) conf = 1;
            if(i>=2&&j<2) conf = 2;
            if(i>=2&&j>=2) conf = 3;
            string dipconf[4] = {"II","IF","FI","FF"};
            
            if(DEBUG){
            cout << "Including "<< dipconf[conf] << " dipole with ";
            cout << PID_to_Name(EXTID[i]) << "(" << i+1 << ") as emitter and ";
            cout << PID_to_Name(EXTID[j]) << "(" << j+1 << ") as spectator" << endl;
//             if(DEBUG>1){
//                 cout << "Emitter momentum = {"<<pe[0]<<","<<pe[1]<<","<<pe[2]<<","<<pe[3]<<","<<pe[4]<<"}"<<endl;
//                 cout << "Spectator momentum = {"<<ps[0]<<","<<ps[1]<<","<<ps[2]<<","<<ps[3]<<","<<ps[4]<<"}"<<endl;
//             }
            }
                        
//             if(conf==0)  II_Poles(pe,ps,ren_scale,aux);
//             if(conf==1) {IF_Poles(pe,ps,ren_scale,aux);ch*=(-1.);}
//             if(conf==2) {FI_Poles(pe,ps,ren_scale,aux);ch*=(-1.);}
//             if(conf==3)  FF_Poles(pe,ps,ren_scale,aux);
            
            if(DEBUG>1) cout << "This dipole contribution = {" << ch*aux[0] << "," << ch*aux[1] << "}" << endl<<endl;
            
//             ret_val[0] += ch*aux[0];
//             ret_val[1] += ch*aux[1];
        
            }
        }
    }
    
    
}








main(int argc, char* argv[])
{
    
//     Process now udbar ---> Wbbar  ///
    
    int next = 5;
    int id[5] = {2,-1,24,5,-5};
    Create_Subtracted(next,id);
    
    
    
    
   
}