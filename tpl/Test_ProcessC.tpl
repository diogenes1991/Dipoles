#include <string>
#include <cstdlib>
#include <ctime>
#include <cmath>

// For C++ programs, the following needs to be included in order to use
// NLOX -  the only include needed to use the generated NLOX processes.
#include "nlox_olp.h"
#include "Phase_Space_Tools.h"

int main(int argc, char* argv[]){
    
    // The first string in argv, argv[0], is the command used to invoke
    // the program, so argc is 5 or 6 for 4 or 5 arguments given.
    if (argc > 6) {
    std::cout << "Too many input arguments!" << std::endl;
    std::cout << "Only five input arguments allowed, i.e." << std::endl;
    std::cout << "  1) the integer that specifies the subprocess: see SUBPROCESSES file," << std::endl;
    std::cout << "  2) the string that specifies the interference type: tree_tree or tree_loop," << std::endl;
    std::cout << "  3) the string that specifies the coupling-power combination: asXaeY or gI'eJ'_gIeJ,"<< std::endl;
    std::cout << "  4) a fixed scale mu in GeV,"<< std::endl;
    std::cout << "  5) an optional 5th argument may be given: an integer 0 or 1 to"<< std::endl;
    std::cout << "  specify whether color-correlated Born MEs should be evaluated"<< std::endl;
    std::cout << "  and printed in addition (1) or not (0; default). Note that for"<< std::endl;
    std::cout << "  tree_loop there are no color-correlated Born MEs generated."<< std::endl;
    std::cout << "Abort now."<< std::endl;
    abort();
    } 
    else if (argc < 5) {
    std::cout << "Not enough input arguments!" << std::endl;
    std::cout << "Four input arguments needed, i.e." << std::endl;
    std::cout << "  1) the integer that specifies the subprocess: see SUBPROCESSES file," << std::endl;
    std::cout << "  2) the string that specifies the interference type: tree_tree or tree_loop," << std::endl;
    std::cout << "  3) the string that specifies the coupling-power combination: asXaeY or gI'eJ'_gIeJ," << std::endl;
    std::cout << "  4) a fixed scale mu in GeV."<< std::endl;
    std::cout << std::endl;
    std::cout << "  An optional 5th argument may be given: an integer 0 or 1 to"<< std::endl;
    std::cout << "  specify whether color-correlated Born MEs should be evaluated"<< std::endl;
    std::cout << "  and printed in addition (1) or not (0; default). Note that for"<< std::endl;
    std::cout << "  tree_loop there are no color-correlated Born MEs generated."<< std::endl;
    std::cout << "Abort now."<< std::endl;
    abort();
    }

    NLOX_OLP_Start(NULL,NULL);
    
    // //////////////////////////////////////////////////////////////////
    // PHASE SPACE POINT
    // (this is absolutely user/MC dependent, and only the form of the PSP
    // pp is relevant)
    // pp stores the phase-space point. It is an array of double precision
    // numbers of length 5 * (number of external particles).
    //
    // The elements of pp are
    //
    //    pp = [p1t, p1x, p1y, p1z, m1, p2t, p2x, ...]
    //
    // This is the format used by the BLHA2 standard. 
    //
    // NOTE: NLOX does currently not use the mass components, nor the last
    // five entries of pp.
    //
    // There may be several phase-space points listed. Comment out all but
    // the one you want to use.
    
    ////phase_space_point
    
    // Convert PSP to BLHA format.
    ////setPCM


    int nextRadi= ####nExtRadi####;
    int nextBorn = nextRadi-1;

    // Center of Mass energy in GeV

    double sqrts = 1000;
    FourVector P(sqrts,0,0,0);
    
    FourVector PIn[2];
    double MIn[2];
    double rIn[2];
    double JIn=1;  

    FourVector PFB[2];
    double MFB[2];
    double rFB[2];
    double JFB=1; 

    FourVector PFR[2];
    double MFR[2];
    double rFR[2];
    double JFR=1;   
        
    ####Evaluate BornsC####

    

    return 0;
}