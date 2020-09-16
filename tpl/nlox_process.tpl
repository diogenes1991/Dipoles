#ifndef __NLOX_PROCESS_H_
#define __NLOX_PROCESS_H_

####Include Born####
####Include Radi####

class Process {
  Subprocess** subproc;
  int numSubprocesses;
  std::unordered_map<std::string,int> AmpMap;

  public:
    ProcessConst pc;
    
    Process() {
      numSubprocesses = ####NSubProcesses####;
      subproc = new Subprocess* [numSubprocesses];

      // Initialize classes derived from Subprocess
        ####Construct Born####
        ####Construct Radi####

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

      *acc = ( evalOK ? ( (IRP[1] == 0.) ? 0 : std::abs((rval[1]/IRP[1])-1. )) : -1. );
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
      
      *acc = ( evalOK ? ( (IRP[1] == 0.) ? 0 : std::abs((rval2[1]/IRP[1])-1. )) : -1. );
    }

    void update_mass(const std::string & name) {
      for (int i = 0; i < numSubprocesses; ++i) {
        subproc[i]->set_mass(name, pc.get_mass(name));
      }
    }

};

#endif // __NLOX_PROCESS_H_
