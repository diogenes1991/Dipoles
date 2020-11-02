#ifndef _PLUS_DISTRIBUTION_H
#define _PLUS_DISTRIBUTION_H

class Plus_Distribution{
    double (*Divergent)(double*), (*Regular)(double*); 
    public:
        
        /// Constructors & De-Constructor ///
        Plus_Distribution();
        Plus_Distribution(double (*)(double*) throw (),double (*)(double*) throw ());
        
        
        /// Functions to set and evaluate the integrands ///
        void set_D(double (*)(double*) throw ());
        void set_R(double (*)(double*) throw ());
        
        double eval_D(double* x);
        double eval_R(double* x);
        
};

#endif
