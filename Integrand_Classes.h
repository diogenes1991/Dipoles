#ifndef _Integrand_Classes_H
#define _Integrand_Classes_H



double ZERO(double x){
    return 0.0;
}

class Plus_Distribution{
    double (*Divergent)(double), (*Regular)(double); 
    public:
        
        /// Constructors ///
        Plus_Distribution();
        Plus_Distribution(double (*)(double) throw (),double (*)(double) throw ());
        void set_D(double (*)(double) throw ());
        void set_R(double (*)(double) throw ());
        
        double eval_D(double x);
        double eval_R(double x);
        
};

Plus_Distribution::Plus_Distribution(){
    Divergent = &ZERO;
    Regular = &ZERO;
}

Plus_Distribution::Plus_Distribution(double (*D)(double) throw (),double (*R)(double) throw ()){
    
    double (*aux1)(double) = D;
    double (*aux2)(double) = R;
    set_D(aux1);
    set_R(aux2);
}

void Plus_Distribution::set_D(double (*D)(double) throw ()) {Divergent = D;}
void Plus_Distribution::set_R(double (*R)(double) throw ()) {Regular = R;}

double Plus_Distribution::eval_D( double x ){
    return Divergent(x);
}
double Plus_Distribution::eval_R( double x ){
    return Regular(x);
}




#endif
