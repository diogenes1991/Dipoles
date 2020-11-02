#ifndef _ARGUMENTS_H
#define _ARGUMENTS_H

struct SubArg{

    std::vector<FourVector> psp;
    std::string cp;
    double* rval;

};

struct PluArg{

    std::string cp;
    std::vector<FourVector> psp;
    double J;
    std::vector<FourVector> psp_1;
    double J_1;
    std::vector<double> mlist;
    double x;
    double* rval

};

struct EndArg{

    std::string cp;
    std::vector<FourVector> psp;
    double* rval

};


#endif
