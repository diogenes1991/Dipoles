
template <class T>
FourMatrixT<T> Boost_Matrix_II(FourVectorT<T> pa, FourVectorT<T> pb, FourVectorT<T> k){
    FourVectorT<T> Pab = pa + pb - k;
    FourVectorT<T> Pab_tilde = pa_II(pa,pb,k) + pb;
    FourVectorT<T> Sum = Pab + Pab_tilde;
    FourMatrixT<T> out;
    FourMatrixT<T> aux = TensorProduct(Pab,Pab_tilde);
                   aux = aux * (2.0/(Pab*Pab));
                   out = out + aux;
                   aux = TensorProduct(Sum,Sum);
                   aux = aux * (1.0/((Sum)*(Pab)));
                   out = out - aux;
    for(int i=0;i<=3;i++) out.M[i][i] += 1.0;
    #if VALIDATE 
        std::cout<<out<<std::endl;
    #endif

    return out;
}
