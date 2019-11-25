#include <iostream>
#include <fstream>
#include <cmath>
#include "Vector_Real.h"
#include "Phase_Spaces_Real.h"
#include "Dipole_Definitions.h"
// #include "IR_Pole_NLOX.h"
#include "qqbar_qqbar.h"
#include "qqbar_qpqpbar.h"
#include <time.h>



main(int argc, char* argv[]){
    
    cout.precision(16);
    
    cout << "Success in loading libraries" <<endl;   
    
    Vector P(1000,0,0,0);
    Vector pi[2];
    double mi[2];
    Vector pf[2];
    double mf[2];
    double ri[2];
    double rf[2];
    double J[2];
    
    for(int i=0;i<=1;i++){
        mi[i]=0;
        mf[i]=0;
        ri[i]=0;
        rf[i]=aleatorio(0,1);
    }
    
    Recursive_PSP(P,2,pi,mi,ri,J[0]);
    Recursive_PSP(P,2,pf,mf,rf,J[0]);
    
    for(int i=0;i<=1;i++){
        cout << "p"<<i+1<<" = "<< pi[i] << endl;
        cout << "p"<<i+3<<" = "<< pf[i] << endl;
    }
    
    double pp[5*4];
    for(int i=0;i<=1;i++){
        
        pp[5*i+0] = pi[i].get_t();
        pp[5*i+1] = pi[i].get_x();
        pp[5*i+2] = pi[i].get_y();
        pp[5*i+3] = pi[i].get_z();
        pp[5*i+4] = mi[i];
        
        pp[5*i+10] = pf[i].get_t();
        pp[5*i+11] = pf[i].get_x();
        pp[5*i+12] = pf[i].get_y();
        pp[5*i+13] = pf[i].get_z();
        pp[5*i+14] = mf[i];
        
    }
    
    
    double ppp[20]={250,0,0,250,0,250,0,0,-250,0,
    250,-19.18429364851453,113.7028340258232,221.8189090490341,0,250,19.18429364851453,-113.7028340258232,-221.8189090490341,0};

    Vector p1(250,0,0,250);
    Vector p2(250,0,0,-250);
    Vector p3(250,-19.18429364851453,113.7028340258232,221.8189090490341);
    Vector p4(250,19.18429364851453,-113.7028340258232,-221.8189090490341);
    
    Vector PPP[4];
    PPP[0]=p1;
    PPP[1]=p2;
    PPP[2]=p3;
    PPP[3]=p4;
    double MMM[4] = {0,0,0,0};
    
    
    
    
    cout << "We proceed to check all the different colour-correlated amplitudes" << endl;
    for(int i=0;i<=3;i++){
        for(int j=i+1;j<=3;j++){
            int col[2] = {i+1,j+1};
            cout << "qqbar_qqbar <T"<<i+1<<"T"<<j+1<<"> = " << qqbar_qqbar_color(PPP,MMM,col)/36.<< endl;
        }
    }
    cout << "Born Matrix Element = " << qqbar_qqbar(PPP,MMM)/36. << endl;
    
    cout << "We proceed to check all the different colour-correlated amplitudes" << endl;
    for(int i=0;i<=3;i++){
        for(int j=i+1;j<=3;j++){
            int col[2] = {i+1,j+1};
            cout << "qqbar_qpqpbar <T"<<i+1<<"T"<<j+1<<"> = " << qqbar_qpqpbar_color(PPP,MMM,col)/36.<< endl;
        }
    }
    cout << "Born Matrix Element = " << qqbar_qpqpbar(PPP,MMM)/36. << endl;
    
    double ren_scale = 500;
    double Single_Pole = 0;
    double Double_Pole = 0;
    for(int i=0;i<=3;i++){
        for(int j=0;j<=3;j++){
          if (i!=j){
            int col[2] = {i+1,j+1};
            Double_Pole += qqbar_qqbar_color(PPP,MMM,col);
            Vector Pi,Pj;
            if(i==0) Pi=p1;
            if(i==1) Pi=p2;
            if(i==2) Pi=p3;
            if(i==3) Pi=p4;
            if(j==0) Pj=p1;
            if(j==1) Pj=p2;
            if(j==2) Pj=p3;
            if(j==3) Pj=p4;          
            Single_Pole -= (3./2. + log(ren_scale*ren_scale/(2*(Pi*Pj))))*qqbar_qqbar_color(PPP,MMM,col);
          }
              
        }
        
    }
    
    
    int PIDS[4]={1,-1,1,-1};
    double rets[2];
//     Get_Process_QCD_Poles_NEW(4,PIDS,ppp,ren_scale,color_born,rets);

    cout << endl << endl;
    cout << "Born Matrix Element = " << qqbar_qqbar(PPP,MMM)/36. << endl;
    cout << "NLOX's Born Matrix Element = " << 269.5672835097302 << endl;
    cout << "Ratio of Borns = " <<  qqbar_qqbar(PPP,MMM)/(36.*269.5672835097302) << endl;
    
    cout << "Double Pole / Born = " << Double_Pole / qqbar_qqbar(PPP,MMM) << endl;
    cout << "NLOX's Double Pole = " << 16/3. << endl;
    cout << "Ratio of DP's = " << Double_Pole/qqbar_qqbar(PPP,MMM) / (16/3.) << endl;
    
    cout << "Single Pole / Born = " << Single_Pole / qqbar_qqbar(PPP,MMM) << endl;
    cout << "NLOX's Single Pole = " << -6.188017708600925 << endl;
    cout << "Ratio of SP's = " << Single_Pole / (-6.188017708600925*qqbar_qqbar(PPP,MMM)) << endl;
    
    
}