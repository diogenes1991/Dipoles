#ifndef _IR_Pole_NLOX_H
#define _IR_Pole_NLOX_H

#include <stdio.h>
#include <stdlib.h>
#include <cmath>

using namespace std;

int DEBUG = 0;

double Lorentz_Dot( double* a, double* b ){
    return a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3];
    
}

double lambda( double x, double y, double z ){
    return x*x + y*y + z*z - 2*x*y - 2*x*z - 2*y*z; 
}

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

int signo_log = 1;

void II_Poles(double *pa, double *pb, double ren_scale, double *ret_val){
    
    double sab = 2*Lorentz_Dot(pa,pb) + pa[4]*pa[4] + pb[4]*pb[4];
    double sabbar = 2*Lorentz_Dot(pa,pb);
    double lambda_ab = lambda(sab,pa[4]*pa[4],pb[4]*pb[4]);
    
    /// Mass-Scenario Decission, only the emitter matters ///     
    if ( pa[4] != 0 ){
    if(DEBUG) cout << "Massive emitter scenario triggered" << endl;
    ret_val[0] = 0;
    ret_val[1] = sabbar + 2*pa[4]*pa[4] - sqrt(lambda_ab);
    ret_val[1] = ret_val[1] / (sabbar + 2*pa[4]*pa[4] + sqrt(lambda_ab));
    if(DEBUG>1){ 
        cout << "sab = " << sab << endl;
        cout << "sabbar = " << sabbar << endl;
        cout << "log_arg = " << ret_val[1] << endl;}        
    ret_val[1] = log(ret_val[1]);
    ret_val[1] = sabbar * ret_val[1] / sqrt(lambda_ab);
    ret_val[1] = 1 + ret_val[1];
    }
    
    else{
    if(DEBUG) cout << "Massless emitter scenario triggered" << endl;
        ret_val[0] = 1;
        ret_val[1] = (ren_scale*ren_scale*sab)/(sabbar*sabbar);
    if(DEBUG>1){ 
        cout << "sab = " << sab << endl;
        cout << "sabbar = " << sabbar << endl;
        cout << "log_arg = " << ret_val[1] << endl;}
    ret_val[1] = signo_log*log(ret_val[1]);
    ret_val[1] = 3./2 + ret_val[1];
        
    }
    
    ret_val[0] = -ret_val[0];
    ret_val[1] = -ret_val[1];
    
    
    
}

void FF_Poles(double* pi, double* pj, double ren_scale, double* ret_val){
    
    double sij = pi[4]*pi[4] + pj[4]*pj[4] + 2*Lorentz_Dot(pi,pj);
    double sijbar = 2*Lorentz_Dot(pi,pj);
    double lambda_ij = lambda(sij,pi[4]*pi[4],pj[4]*pj[4]);
    
    ///Mass Scenarios ///
    if ( pi[4] != 0 ){
    if(DEBUG) cout << "Massive emitter scenario triggered" << endl;
    ret_val[0] = 0;
    ret_val[1] = sijbar + 2*pi[4]*pi[4] - sqrt(lambda_ij);
    ret_val[1] = ret_val[1] / (sijbar + 2*pi[4]*pi[4] + sqrt(lambda_ij));
    if(DEBUG>1){
        cout << "sij = " << sij << endl;
        cout << "sijbar = " << sijbar << endl;
        cout << "log_arg = " << ret_val[1] << endl;}
    ret_val[1] = log(ret_val[1]);
    ret_val[1] = sijbar * ret_val[1] / sqrt(lambda_ij);
    ret_val[1] = 1 + ret_val[1];
    }
    
    else{
    if(DEBUG) cout << "Massless emitter scenario triggered" << endl;
    ret_val[0] = 1;
    ret_val[1] = (ren_scale*ren_scale*sij)/(sijbar*sijbar);
    if(DEBUG>1){
        cout << "sij = " << sij << endl;
        cout << "sijbar = " << sijbar << endl;
        cout << "log_arg = " << ret_val[1] << endl;}
    ret_val[1] = signo_log*log(ret_val[1]);
    ret_val[1] = 3./2 + ret_val[1];
    
        
    }
    
    ret_val[0] = -ret_val[0];
    ret_val[1] = -ret_val[1];
    
}

void IF_Poles(double* pa, double* pi, double ren_scale, double* ret_val){
    double tai = pa[4]*pa[4] + pi[4]*pi[4] - 2*Lorentz_Dot(pa,pi);
    double taibar = -2*Lorentz_Dot(pa,pi);
    double lambda_ai = lambda(tai,pa[4]*pa[4],pi[4]*pi[4]);
    
     /// Mass-Scenario Decission, only the emitter matters ///     
    if ( pa[4] != 0 ){
    if(DEBUG) cout << "Massive emitter scenario triggered" << endl;
    ret_val[0] = 0;
    ret_val[1] = taibar - 2*pa[4]*pa[4] - sqrt(lambda_ai);
    ret_val[1] = ret_val[1] / (taibar - 2*pa[4]*pa[4] + sqrt(lambda_ai));
    if(DEBUG>1){
        cout << "tai = " << tai << endl;
        cout << "taibar = " << taibar << endl;
        cout << "log_arg = " << ret_val[1] << endl;}
    ret_val[1] = log(ret_val[1]);
    ret_val[1] = taibar * ret_val[1] / sqrt(lambda_ai);
    ret_val[1] = 1 + ret_val[1];
    }
    
    else{
    if(DEBUG) cout << "Massless emitter scenario triggered" << endl;
    ret_val[0] = 1;
    ret_val[1] = -(ren_scale*ren_scale*tai)/(taibar*taibar);
    if(DEBUG>1){ 
        cout << "tai = " << tai << endl;
        cout << "taibar = " << taibar << endl;
        cout << "log_arg = " << ret_val[1] << endl;}
    ret_val[1] = signo_log*log(ret_val[1]);
    ret_val[1] = 3./2 + ret_val[1];
        
    }
   
   ret_val[0] = -ret_val[0];
   ret_val[1] = -ret_val[1];
    
    
}

void FI_Poles(double* pi, double* pa, double ren_scale, double* ret_val){
    double tia = pa[4]*pa[4] + pi[4]*pi[4] - 2*Lorentz_Dot(pa,pi);
    double tiabar = -2*Lorentz_Dot(pa,pi);
    double lambda_ia = lambda(tia,pa[4]*pa[4],pi[4]*pi[4]);
    
     /// Mass-Scenario Decission, only the emitter matters ///     
    if ( pi[4] != 0 ){
    if(DEBUG) cout << "Massive emitter scenario triggered" << endl;
    ret_val[0] = 0;
    ret_val[1] = tiabar - 2*pi[4]*pi[4] - sqrt(lambda_ia);
    ret_val[1] = ret_val[1] / (tiabar - 2*pi[4]*pi[4] + sqrt(lambda_ia));
    if(DEBUG>1){
        cout << "tia = " << tia << endl;
        cout << "tiabar = " << tiabar << endl;
        cout << "log_arg = " << ret_val[1] << endl;}
    ret_val[1] = log(ret_val[1]);
    ret_val[1] = tiabar * ret_val[1] / sqrt(lambda_ia);
    ret_val[1] = 1 + ret_val[1];
    }
    
    else{
    if(DEBUG) cout << "Massless emitter scenario triggered" << endl;
        ret_val[0] = 1;
        ret_val[1] = -(tia*ren_scale*ren_scale)/(tiabar*tiabar);
//         ret_val[1] = -(ren_scale*ren_scale)/(tiabar);
    if(DEBUG>1){ 
        cout << "tia = " << tia << endl;
        cout << "tiabar = " << tiabar << endl;
        cout << "log_arg = " << ret_val[1] << endl;}
//         ret_val[1] = signo_log*log(ret_val[1]);
        ret_val[1] = log(ret_val[1]);
        ret_val[1] = 3./2 + ret_val[1];
        
    }
   
   ret_val[0] = -ret_val[0];
   ret_val[1] = -ret_val[1];
    
    
}

void Get_Process_EW_Poles( const int next, int* subproc_id, double *pp, double ren_scale, double* ret_val){
    /// Returns the Pole for a defined process ///

    if(DEBUG){
        for(int i=0;i<=next-1;i++){
    if(i==0) cout << "QED dipole routine triggered for the process " << endl;
    if(i==2) cout << " ---> ";
    cout << PID_to_Name(subproc_id[i]) << "(" << i+1 << ")";
    if(i!=next-1&&i!=1) cout << " + ";
    else if(i!=1) cout << endl<<endl;
    }
    }

    double aux[2];
    
    bool USEALTCH = 0;
    
    for(int i=0;i<=next-1;i++){
        for(int j=0;j<=next-1;j++){
            double ch;
            if(USEALTCH) ch = Charge(subproc_id[i])*Charge(subproc_id[j])*signo(subproc_id[i])*signo(subproc_id[j]);
            else ch = Charge(abs(subproc_id[i]))*Charge(abs(subproc_id[j]))*signo(subproc_id[i])*signo(subproc_id[j]);
            if(ch&&j!=i){
            double pe[5],ps[5];
            for(int k=0;k<5;k++){pe[k]=pp[5*i+k];ps[k]=pp[5*j+k];}
            int conf = 0;
            if(i<2&&j<2) conf = 0;
            if(i<2&&j>=2) conf = 1;
            if(i>=2&&j<2) conf = 2;
            if(i>=2&&j>=2) conf = 3;
            string dipconf[4] = {"II","IF","FI","FF"};
            
            if(DEBUG){
            cout << "Including "<< dipconf[conf] << " dipole with ";
            cout << PID_to_Name(subproc_id[i]) << "(" << i+1 << ") as emitter and ";
            cout << PID_to_Name(subproc_id[j]) << "(" << j+1 << ") as spectator" << endl;
            if(DEBUG>1){
                cout << "Emitter momentum = {"<<pe[0]<<","<<pe[1]<<","<<pe[2]<<","<<pe[3]<<","<<pe[4]<<"}"<<endl;
                cout << "Spectator momentum = {"<<ps[0]<<","<<ps[1]<<","<<ps[2]<<","<<ps[3]<<","<<ps[4]<<"}"<<endl;
            }
            }
                        
            if(conf==0)  II_Poles(pe,ps,ren_scale,aux);
            if(conf==1) {IF_Poles(pe,ps,ren_scale,aux);ch*=(-1.);}
            if(conf==2) {FI_Poles(pe,ps,ren_scale,aux);ch*=(-1.);}
            if(conf==3)  FF_Poles(pe,ps,ren_scale,aux);
            
            if(DEBUG>1) cout << "This dipole contribution = {" << ch*aux[0] << "," << ch*aux[1] << "}" << endl<<endl;
            
            ret_val[0] += ch*aux[0];
            ret_val[1] += ch*aux[1];
        
            }
        }
    }
    
}

void IFFI_Poles(double* pa, double* pi, double ren_scale, double* ret_val){
    double tia = pa[4]*pa[4] + pi[4]*pi[4] - 2*Lorentz_Dot(pa,pi);
    double tiabar = -2*Lorentz_Dot(pa,pi);
    double lambda_ia = lambda(tia,pa[4]*pa[4],pi[4]*pi[4]);
    
     /// Mass-Scenario Decission configuration detection ///
    
    int config;
    
    /// 0 --> mi=0 ma=0      ///
    /// 1 --> mi=0 ma!=0     ///
    /// 2 --> mi!=0 ma=0     ///
    /// 3 --> mi!=0 ma!=0    ///
    
    if(pi[4]==0){
        if(pa[4]==0) config=0;
        else config=1;}
    else{
        if(pa[4]==0) config=2;
        else config=3;}
        
        
    if (config==3){
        if(DEBUG) cout << "Fully massive scenario triggered" << endl;
        
        ret_val[0] = 0;
        if(DEBUG>1){ 
            cout << "tia = " << tia << endl;
            cout << "tiabar = " << tiabar << endl;
            cout << "log_arg = " << ret_val[1] << endl;}
        ret_val[1] = log((tiabar - 2*(pi[4]*pi[4]) - sqrt(lambda_ia))/(tiabar - 2*(pi[4]*pi[4]) + sqrt(lambda_ia)));
        ret_val[1] = ret_val[1] + log((tiabar - 2*(pa[4]*pa[4]) - sqrt(lambda_ia))/(tiabar - 2*(pa[4]*pa[4]) + sqrt(lambda_ia)));
        ret_val[1] = tiabar * ret_val[1] / sqrt(lambda_ia);
        ret_val[1] = 2 + ret_val[1];        
    }
    
    
    if (config==0){
        if(DEBUG) cout << "Fully massless scenario triggered" << endl;
        ret_val[0] = 2;
        ret_val[1] = -ren_scale*ren_scale;
        ret_val[1] = ret_val[1]/tiabar;
        ret_val[1] = 2*log(ret_val[1]);
        ret_val[1] = 3 + ret_val[1];
    }
    
    
    
    if (config==1){
        if(DEBUG) cout << "Massless final-State massive initial-state scenario triggered" << endl;
        ret_val[0] = 1;
        ret_val[1] = -ren_scale*pa[4];
        ret_val[1] = ret_val[1]/tiabar;
        ret_val[1] = 2*log(ret_val[1]);
        ret_val[1] = 5/2. + ret_val[1];
    }
    
    
    if(config==2){
        if(DEBUG) cout << "Massless initial-State massive final-state scenario triggered" << endl;
        ret_val[0] = 1;
        ret_val[1] = ren_scale*pi[4];
        ret_val[1] = ret_val[1]/(-tiabar);
        ret_val[1] = 2*log(ret_val[1]);
        ret_val[1] = 5/2. + ret_val[1];
    }
        
   ret_val[0] = -ret_val[0];
   ret_val[1] = -ret_val[1];
    
    
}

void Get_Process_EW_Poles_NEW( const int next, int* subproc_id, double *pp, double ren_scale, double* ccborn, double* ret_val){
    
    int Scheme=1;
    ret_val[0]=0;
    ret_val[1]=0;
    
    if(Scheme!=0&&Scheme!=1){
        cout << "Error: Scheme specification not detected "<<endl;
        abort();
    }
    /// Returns the Pole for a defined process ///

    if(DEBUG){
        for(int i=0;i<=next-1;i++){
    if(i==0) cout << "NEW QED dipole routine triggered for the process " << endl;
    if(i==2) cout << " ---> ";
    cout << PID_to_Name(subproc_id[i]) << "(" << i+1 << ")";
    if(i!=next-1&&i!=1) cout << " + ";
    else if(i!=1) cout << endl<<endl;
    }
    if(Scheme==1) cout << "Gmu Scheme active" << endl;
    if(Scheme==0) cout << "Alpha0 Scheme active, all fermion masses should be non-zero!" << endl;
    }

    double aux[2];
    
    bool USEALTCH = 0;
    
    /// f,W --> f+A , W+A Routine ///
    
    for(int i=0;i<=next-1;i++){
        for(int j=0;j<=next-1;j++){
            double ch;
            if(USEALTCH) ch = Charge(subproc_id[i])*Charge(subproc_id[j])*signo(subproc_id[i])*signo(subproc_id[j]);
            else ch = Charge(abs(subproc_id[i]))*Charge(abs(subproc_id[j]))*signo(subproc_id[i])*signo(subproc_id[j]);
            if(ch&&j!=i){
            double pe[5],ps[5];
            for(int k=0;k<5;k++){pe[k]=pp[5*i+k];ps[k]=pp[5*j+k];}
            int conf = 0;
            if(i<2&&j<2) conf = 0;
            if(i<2&&j>=2) conf = 1;
            if(i>=2&&j<2) conf = 2;
            if(i>=2&&j>=2) conf = 3;
            string dipconf[4] = {"II","IF","FI","FF"};
            
            if(DEBUG){
            cout << "Including "<< dipconf[conf] << " dipole with ";
            cout << PID_to_Name(subproc_id[i]) << "(" << i+1 << ") as emitter and ";
            cout << PID_to_Name(subproc_id[j]) << "(" << j+1 << ") as spectator" << endl;
            if(DEBUG>1){
                cout << "Emitter momentum = {"<<pe[0]<<","<<pe[1]<<","<<pe[2]<<","<<pe[3]<<","<<pe[4]<<"}"<<endl;
                cout << "Spectator momentum = {"<<ps[0]<<","<<ps[1]<<","<<ps[2]<<","<<ps[3]<<","<<ps[4]<<"}"<<endl;
            }
            }
                        
            if(conf==0)  II_Poles(pe,ps,ren_scale,aux);
            if(conf==1) {IFFI_Poles(pe,ps,ren_scale,aux);ch*=(-1./2);}
            if(conf==2) {IFFI_Poles(pe,ps,ren_scale,aux);ch*=(-1./2);}
            if(conf==3)  FF_Poles(pe,ps,ren_scale,aux);
            
            if(DEBUG>1) cout << "This dipole contribution = {" << ch*aux[0] << "," << ch*aux[1] << "}" << endl<<endl;
            
            ret_val[0] += ch*aux[0];
            ret_val[1] += ch*aux[1];
        
            }
        }
    }
    
    /// A --> ffbar + PDF Counterterm ///
    if(Scheme==1){
    for(int i=0;i<=next-1;i++){
        if(subproc_id[i]==22)ret_val[1] -= 40./9;
        }
    }    
   ret_val[0] *=ccborn[0];
   ret_val[1] *=ccborn[0];
    
    
}

void Get_Process_QCD_Poles( const int next, int* subproc_id, double* pp, double ren_scale,double *born, double* ret_val){
    /// Returns the Pole for a defined process ///
    double CF = 4./3;

    if(DEBUG){
        for(int i=0;i<=next-1;i++){
    if(i==0) cout << "QCD dipole routine triggered for the process " << endl;
    if(i==2) cout << " ---> ";
    cout << PID_to_Name(subproc_id[i]) << "(" << i+1 << ")";
    if(i!=next-1&&i!=1) cout << " + ";
    else if(i!=1) cout << endl<<endl;
    }
    }
    
    double aux[2];
    
    for(int i=0;i<=next-1;i++){
        for(int j=0;j<=next-1;j++){
            double ch = CF*IsColored(subproc_id[i])*IsColored(subproc_id[j])*signo(subproc_id[i])*signo(subproc_id[j]);
            if(ch&&j!=i){
            double pe[5],ps[5];
            for(int k=0;k<5;k++){pe[k]=pp[5*i+k];ps[k]=pp[5*j+k];}
            int conf = 0;
            if(i<2&&j<2) conf = 0;
            if(i<2&&j>=2) conf = 1;
            if(i>=2&&j<2) conf = 2;
            if(i>=2&&j>=2) conf = 3;
            string dipconf[4] = {"II","IF","FI","FF"};
            
            if(DEBUG){
            cout << "Including "<< dipconf[conf] << " dipole with ";
            cout << PID_to_Name(subproc_id[i]) << "(" << i+1 << ") as emitter and ";
            cout << PID_to_Name(subproc_id[j]) << "(" << j+1 << ") as spectator" << endl;
            if(DEBUG>2){
                cout << "Emitter momentum = {"<<pe[0]<<","<<pe[1]<<","<<pe[2]<<","<<pe[3]<<","<<pe[4]<<"}"<<endl;
                cout << "Spectator momentum = {"<<ps[0]<<","<<ps[1]<<","<<ps[2]<<","<<ps[3]<<","<<ps[4]<<"}"<<endl;
            }
            }
                        
            if(conf==0) II_Poles(pe,ps,ren_scale,aux);
            if(conf==1) {IF_Poles(pe,ps,ren_scale,aux);ch*=(-1);}
            if(conf==2) {FI_Poles(pe,ps,ren_scale,aux);;ch*=(-1);}
            if(conf==3) FF_Poles(pe,ps,ren_scale,aux);
            
            if(DEBUG>1) cout << "This dipole contribution = {" << ch*aux[0] << "," << ch*aux[1] << "}" << endl<<endl;
            
            ret_val[0] += ch*aux[0];
            ret_val[1] += ch*aux[1];
        
            }
        }
    }
    
}


double gamma_I(int PID){
    if(IsColored(PID)){
        if(PID!=21) return 3./2;
        else return 11./6 - 5./9;
    }
    else return 0;
    
}

void Get_Process_QCD_Poles_NEW( const int next, int* subproc_id, double* pp, double ren_scale, double *ccborn, double* ret_val){

    ret_val[0] = 0;
    ret_val[1] = 0;
    if(DEBUG){
        for(int i=0;i<=next-1;i++){
    if(i==0) cout << "NEW QCD dipole routine triggered for the process " << endl;
    if(i==2) cout << " ---> ";
    cout << PID_to_Name(subproc_id[i]) << "(" << i+1 << ")";
    if(i!=next-1&&i!=1) cout << " + ";
    else if(i!=1) cout << endl;
    }
    }
    
    
    
    int CC_count = 1;
    for(int i=0;i<=next-1;i++){
        for(int j=i+1;j<=next-1;j++){
            if(IsColored(subproc_id[i])&&IsColored(subproc_id[j])){
            double pe[5],ps[5];
            for(int k=0;k<5;k++){pe[k]=pp[5*i+k];ps[k]=pp[5*j+k];}            
            if(DEBUG){
                cout << "Now including the <T"<<i+1<<"T"<<j+1<<"> = "<< ccborn[CC_count] << endl;
            if(DEBUG>1){
                cout << "Emitter momentum = {"<<pe[0]<<","<<pe[1]<<","<<pe[2]<<","<<pe[3]<<","<<pe[4]<<"}"<<endl;
                cout << "Spectator momentum = {"<<ps[0]<<","<<ps[1]<<","<<ps[2]<<","<<ps[3]<<","<<ps[4]<<"}"<<endl;
                }
            }
                        
            

            ret_val[0] -= 2*ccborn[CC_count];
            ret_val[1] -= (gamma_I(subproc_id[i]) + gamma_I(subproc_id[j]) + 2*log(ren_scale*ren_scale/(2*(Lorentz_Dot(pe,ps)))))*ccborn[CC_count];
            
            }
            CC_count += 1;
        }
    }
    
//     for(int i=0;i<=next-1;i++){
//         if(subproc_id[i]==21)ret_val[1] -= 5./3*ccborn[0];
//         }
    
}

#endif
