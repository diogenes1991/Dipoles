#ifndef __NLOX_PROCESS_H__
#define __NLOX_PROCESS_H__

#include "Phase_Space_Tools.h"

####Include Channels####

class Process {
  Subprocess** subproc;
  int numSubprocesses;

  public:
    ProcessConst pc;
    std::unordered_map<std::string,int> AmpMap;
    
    Process() {
      numSubprocesses = ####NSubProcesses####;
      subproc = new Subprocess* [numSubprocesses];

      // Initialize classes derived from Subprocess
        ####Construct Channels####

        ####Amplitude Map####
        
    }

    ~Process() {
      for (int i = 0; i != numSubprocesses; i++)
        delete subproc[i];
      delete [] subproc;
    }
    
    void evaluate_all(int i, double* pp, int next, double mu, double* rval, double* acc) {
      const double virtfac = -1., CTfac = 1.,IRfac = 1.;
      double virtc2, virtc1, virtc0;
      double CTc2, CTc1, CTc0;
      double IRP[2]={0.,0.};
      
      bool evalOK = true;

      // const double nlofac = 1.;
      // const double PI  = 3.14159265358979323846264338327950;
      const double nlofac = 1./8./M_PI/M_PI;

      std::vector<FourVector> psp;
      for ( int l = 0; l < (next - 1)*5; l+=5 ) {
        FourVector p(pp[l],pp[l+1],pp[l+2],pp[l+3]);
        psp.push_back(p);
      }

      subproc[i]->setPSP2(psp);

      evalOK = subproc[i]->virt(mu, &virtc2, &virtc1, &virtc0, virtfac);
      subproc[i]->CT(mu, &CTc2, &CTc1, &CTc0, CTfac);
      subproc[i]->poles(mu,IRP,IRfac);
      
      rval[0] = nlofac*(virtc2 + CTc2); // Double pole
      rval[1] = nlofac*(virtc1 + CTc1); // Single pole
      rval[2] = nlofac*(virtc0 + CTc0); // Finite
      rval[3] = subproc[i]->born();     // |Born|^2

      *acc = ( evalOK ? ( (IRP[1] == 0.) ? 0 : std::abs((rval[1]/IRP[1])+1. )) : -1. );
    }

    void evaluate(int i, std::string type, std::string cp, double* pp, int next, double mu, double* rval2, double* acc) {
      const double virtfac = -1., CTfac = 1.;
      double virtc2, virtc1, virtc0;
      double CTc2, CTc1, CTc0;
      
      bool evalOK = true;

      // const double nlofac = 1.;
      // const double PI  = 3.14159265358979323846264338327950;
      const double nlofac = 1./8./M_PI/M_PI;

      std::vector<FourVector> psp;
      for ( int l = 0; l < (next - 1)*5; l+=5 ) {
        FourVector p(pp[l],pp[l+1],pp[l+2],pp[l+3]);
        psp.push_back(p);
      }

      subproc[i]->setPSP2(psp);

      if ( type == "tree_tree" ) {
        rval2[0] = 0.;                    // tree_tree double pole
        rval2[1] = 0.;                    // tree_tree single pole
        rval2[2] = subproc[i]->born2(cp); // tree_tree finite
      } else if ( type == "tree_loop" ) {
        evalOK = subproc[i]->virt2(cp, mu, &virtc2, &virtc1, &virtc0, virtfac);
        subproc[i]->CT2(cp, mu, &CTc2, &CTc1, &CTc0, CTfac);
        rval2[0] = nlofac*(virtc2 + CTc2); // tree_loop double pole
        rval2[1] = nlofac*(virtc1 + CTc1); // tree_loop single pole
        rval2[2] = nlofac*(virtc0 + CTc0); // tree_loop finite
      } else {
        std::cout << "Error:" << std::endl;
        std::cout << "The specified intereference type \"" << type << "\" has probably been mistyped for subprocess " << i << "." << std::endl;
        std::cout << "The possible interference types are \"tree_tree\" or \"tree_loop\". Please check for typos." << std::endl;
        std::cout << "Please also check for availability of the desired interference-type for subprocess " << i << "." << std::endl;
        std::cout << "If not available the returned result will be zero." << std::endl;
        abort();
      }
      
      *acc = ( evalOK ? 0 : -1. );
    }
    
    void evaluate_alpha(int i, std::string type, std::string cp, double* pp, int next, double mu, double* rval2, double* acc) {
      const double virtfac = -1., CTfac = 1., IRfac = 1.;
      double virtc2, virtc1, virtc0;
      double CTc2, CTc1, CTc0;
      double IRP[2]={0.,0.};
      
      bool evalOK = true;

      // const double nlofac = 1.;
      // const double PI  = 3.14159265358979323846264338327950;
      const double nlofac = 1./8./M_PI/M_PI;

      std::vector<FourVector> psp;
      for ( int l = 0; l < (next - 1)*5; l+=5 ) {
        FourVector p(pp[l],pp[l+1],pp[l+2],pp[l+3]);
        psp.push_back(p);
      }
      

      subproc[i]->setPSP2(psp);

      if ( type == "tree_tree" ) {
        rval2[0] = 0.;                    // tree_tree double pole
        rval2[1] = 0.;                    // tree_tree single pole
        rval2[2] = subproc[i]->born2Alpha(cp); // tree_tree finite
      } 
      else if ( type == "tree_loop" ){
        evalOK = subproc[i]->virt2Alpha(cp, mu, &virtc2, &virtc1, &virtc0, virtfac);
        subproc[i]->CT2Alpha(cp, mu, &CTc2, &CTc1, &CTc0, CTfac);
        subproc[i]->poles2Alpha(cp,mu,IRP,IRfac);
        rval2[0] = nlofac*(virtc2 + CTc2); // tree_loop double pole
        rval2[1] = nlofac*(virtc1 + CTc1); // tree_loop single pole
        rval2[2] = nlofac*(virtc0 + CTc0); // tree_loop finite
      }
        
      else {
        std::cout << "Error:" << std::endl;
        std::cout << "The specified intereference type \"" << type << "\" has probably been mistyped for subprocess " << i << "." << std::endl;
        std::cout << "The possible interference types are \"tree_tree\" or \"tree_loop\". Please check for typos." << std::endl;
        std::cout << "Please also check for availability of the desired interference-type for subprocess " << i << "." << std::endl;
        std::cout << "If not available the returned result will be zero." << std::endl;
        abort();
      }
      
      *acc = ( evalOK ? ( (IRP[1] == 0.) ? 0 : std::abs((rval2[1]/IRP[1])+1. )) : -1. );
    }

    void evaluate_alpha_sc(int i, std::string type, std::string cp, double* pp, int next, int pn, double** Bmunu){

      std::vector<FourVector> psp;
      for ( int l = 0; l < (next - 1)*5; l+=5 ) {
        FourVector p(pp[l],pp[l+1],pp[l+2],pp[l+3]);
        psp.push_back(p);
      }
      
      subproc[i]->setPSP2(psp);
      
      double cost = pp[5*pn+3]/pp[5*pn];
      double cosp,sinp;
      if(1-std::abs(cost)<1E-17) {cosp = 1;sinp=0.;}
      else{ 
          cosp = pp[5*pn+1]/pp[5*pn]/sqrt(1.-cost*cost);
          sinp = pp[5*pn+2]/pp[5*pn]/sqrt(1.-cost*cost);
      }
      
      FourVector Pth(0,cost*cosp,cost*sinp,-sqrt(1.-cost*cost));
      FourVector Pph(0,-sinp,cosp,0);
      
      subproc[i]->setPol(pn,Pth,Pph);
      double MMCos = subproc[i]->born2Alpha(cp);
      subproc[i]->resetPol();
      
      subproc[i]->setPol(pn,Pth,Pth);
      double MthSq = subproc[i]->born2Alpha(cp);
      subproc[i]->resetPol();
      
      subproc[i]->setPol(pn,Pph,Pph);
      double MphSq = subproc[i]->born2Alpha(cp);
      subproc[i]->resetPol();

      FourMatrixT<double> BmunuM;
      BmunuM = MthSq*TensorProduct(Pth,Pth) + MphSq*TensorProduct(Pph,Pph) + MMCos*(TensorProduct(Pth,Pph)+TensorProduct(Pph,Pth));

      #if VERB
        std::cout<<"|Mph||Mth|cos(a) = "<<MMCos<<std::endl; 
        std::cout<<"      |Mph|^2    = "<<MphSq<<std::endl; 
        std::cout<<"      |Mth|^2    = "<<MthSq<<std::endl; 
        std::cout<<"B_munu = \n"<<BmunuM<<std::endl;
      #endif
      
      
      for(int j=0;j<4;j++){
        for(int k=0;k<4;k++){
          Bmunu[j][k] = BmunuM.M[j][k];
        }
      }

    }

    void evaluate_alpha_cc(int i, std::string type, std::string cp, double* pp, int next, double* Bij){

      std::vector<FourVector> psp;
      for ( int l = 0; l < (next - 1)*5; l+=5 ) {
        FourVector p(pp[l],pp[l+1],pp[l+2],pp[l+3]);
        psp.push_back(p);
      }
      subproc[i]->setPSP2(psp);

      if ( type == "tree_tree" ) {
        subproc[i]->born2Alpha_cc(cp,Bij);
      } 

      else{
        std::cout << "Error:" << std::endl;
        std::cout << "The specified intereference type \"" << type << "\" has probably been mistyped for subprocess " << i << "." << std::endl;
        std::cout << "For Color-Correlated amplitudes only \"tree_tree\" interference types are built." <<std::endl;
        abort();
      }

    }

    void evaluate_alpha_cc_and_sc(int i, std::string type, std::string cp, double* pp, int next, int pn, double*** Bmunuij){

      std::vector<FourVector> psp;
      for ( int l = 0; l < (next - 1)*5; l+=5 ) {
        FourVector p(pp[l],pp[l+1],pp[l+2],pp[l+3]);
        psp.push_back(p);
      }
      
      subproc[i]->setPSP2(psp);

      int NCC = 1+next*(next-1)/2;

      double Bijtt[NCC],Bijpp[NCC],Bijpt[NCC]; 

      double cost = pp[5*pn+3]/pp[5*pn];
      double cosp,sinp;
      if(1-std::abs(cost)<1E-17) {cosp = 1;sinp=0.;}
      else{ 
          cosp = pp[5*pn+1]/pp[5*pn]/sqrt(1.-cost*cost);
          sinp = pp[5*pn+2]/pp[5*pn]/sqrt(1.-cost*cost);
      }
      
      FourVector Pth(0,cost*cosp,cost*sinp,-sqrt(1.-cost*cost));
      FourVector Pph(0,-sinp,cosp,0);
      
      subproc[i]->setPol(pn,Pth,Pph);
      subproc[i]->born2Alpha_cc(cp,Bijpt);
      subproc[i]->resetPol();
      
      subproc[i]->setPol(pn,Pth,Pth);
      subproc[i]->born2Alpha_cc(cp,Bijtt);
      subproc[i]->resetPol();
      
      subproc[i]->setPol(pn,Pph,Pph);
      subproc[i]->born2Alpha_cc(cp,Bijpp);
      subproc[i]->resetPol();
      
      FourMatrixT<double> Utt = TensorProduct(Pth,Pth);
      FourMatrixT<double> Upp = TensorProduct(Pph,Pph);
      FourMatrixT<double> Upt = TensorProduct(Pph,Pth) + TensorProduct(Pth,Pph);

      
      for (int cc=0;cc<NCC;cc++){
        FourMatrixT<double> BmunuijM = Bijtt[cc]*Utt + Bijpp[cc]*Upp + Bijpt[cc]*Upt;
        #if VERB
        std::cout << "cccounter = " << cc << std::endl;
        std::cout<<"|Mph||Mth|cos(a) = "<<Bijpt[cc]<<std::endl; 
        std::cout<<"      |Mph|^2    = "<<Bijpp[cc]<<std::endl; 
        std::cout<<"      |Mth|^2    = "<<Bijtt[cc]<<std::endl; 
        std::cout << BmunuijM << std::endl;
        #endif
        for(int mu=0;mu<4;mu++){
          for(int nu=0;nu<4;nu++){
            Bmunuij[cc][mu][nu] = BmunuijM.M[mu][nu];
          }
        }
      }

    }
    
    void update_mass(const std::string & name) {
      for (int i = 0; i < numSubprocesses; ++i) {
        subproc[i]->set_mass(name, pc.get_mass(name));
      }
    }

};

#endif // __NLOX_PROCESS_H_
