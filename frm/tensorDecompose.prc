#procedure tensorDecompose()
********************************************************************************
*
* tensorDecompose
*
* Tensor decomposition of KK(k(v1))*KK(k(v2))*...*TI(pts, ...)
* (implemented up to rank = 4)
*
********************************************************************************

b KK, TI;
.sort

#message >> tensor decomposition
repeat id KK(k(v1?))*TI(n1?, ?args) = TI(n1, k(v1),?args);
repeat id TI(n?, ?args, k(v1?)) = TI(n,k(v1),?args);
.sort

i c1,c2,c3,c4,c5;
cf f, ff, gg, pp, fi, fj, ppp, gggg, ddd(symmetric);
cf sf;
Keep Brackets;


* rank 1
id TI(n1?,k(v1?), p1?!{k}, ?args) = 
    sum_(c1,1,n1-1,f(c1, k(v1), p1, ?args)*TI(n1,c1, p1, ?args));

repeat id once f(n1?{,>1}, k(v1?), p1?!{k}, p2?, ?args) = f(n1-1, k(v1), p1+p2, ?args);
id f(1,k(v1?), p1?!{k}, ?args) = p1(v1);
b TI;
.sort
Keep Brackets;

* rank 2
id TI(n1?,k(v1?),k(v2?), p1?!{k}, ?args) = 
    sum_(c1,1,n1-1, 
        sum_(c2,1,n1-1,
             f(c1, k(v1), p1, ?args)
            *f(c2, k(v2), p1, ?args)
            *TI(n1,sf(c1,c2), p1, ?args)
        )
    )  
    + d_(v1,v2)*TI(n1,0,0,p1,?args);


repeat id once f(n1?{,>1}, k(v1?), p1?!{k}, p2?, ?args) = f(n1-1, k(v1), p1 + p2, ?args);
id f(1,k(v1?), p1?!{k}, ?args) = p1(v1);

b TI;
.sort
Keep Brackets;

* rank 3
id TI(n1?,k(v1?),k(v2?), k(v3?), p1?!{k}, ?args) = 
    sum_(c1,1,n1-1, 
        sum_(c2,1,n1-1,
            sum_(c3,1,n1-1,
                 f(c1, k(v1), p1, ?args)
                *f(c2, k(v2), p1, ?args)
                *f(c3, k(v3), p1, ?args)
                *TI(n1,sf(c1,c2,c3), p1, ?args)
            )
        )
    ) 
    + sum_(c1,1,n1-1, 
        ff(ddd(v1,v2),f(c1, k(v3), p1, ?args))*TI(n1,0,0,c1,p1,?args)
    );

repeat id once f(n1?{,>1}, k(v1?), p1?!{k}, p2?, ?args) = f(n1-1, k(v1), p1+p2, ?args);
id f(1,k(v1?), p1?!{k}, ?args) = p1(v1);

argument;
    repeat id once f(n1?{,>1}, k(v1?), p1?!{k}, p2?, ?args) = f(n1-1, k(v1), p1 + p2, ?args);
    id f(1,k(v1?), p1?!{k}, ?args) = f(p1,v1);
endargument;

id ff(ddd(v1?,v2?),f(p1?,v3?)) = ff(gg,pp)*fi(v1,v2,v3)*fj(p1);
id ff(gg,pp)*fi(?args2) = distrib_(+1,2,gg,pp,?args2);
id pp(v1?)*fj(p1?) = p1(v1);
id gg(v1?,v2?) = d_(v1,v2);

b TI;
.sort
Keep Brackets;
* rank 4
id TI(n1?,k(v1?),k(v2?), k(v3?), k(v4?), p1?!{k}, ?args) = 
    sum_(c1,1,n1-1, 
        sum_(c2,1,n1-1,
            sum_(c3,1,n1-1,
                sum_(c4,1,n1-1,
                     f(c1, k(v1), p1, ?args)
                    *f(c2, k(v2), p1, ?args)
                    *f(c3, k(v3), p1, ?args)
                    *f(c4, k(v4), p1, ?args)
                    *TI(n1,sf(c1,c2,c3,c4), p1, ?args)
                )
            )
        )
    ) + sum_(c1, 1, n1-1,
            sum_(c2,1,n1-1,
                 ff(ddd(v1,v2),f(c1,k(v3), p1, ?args),f(c2,k(v4), p1, ?args))
                *TI(n1,0,0,sf(c1,c2),p1,?args)
            )
        ) 
    + ff(gg,gg)*fi(v1,v2,v3,v4)*TI(n1,0,0,0,0,p1,?args);

repeat id once f(n1?{,>1}, k(v1?), p1?!{k}, p2?, ?args) = f(n1-1, k(v1), p1+p2, ?args);
id f(1,k(v1?), p1?!{k}, ?args) = p1(v1);

argument;
    repeat id once f(n1?{,>1}, k(v1?), p1?!{k}, p2?, ?args) = f(n1-1, k(v1), p1+p2, ?args);
    id f(1,k(v1?), p1?!{k}, ?args) = f(p1,v1);
endargument;

id ff(ddd(v1?,v2?),f(p1?,v3?), f(p2?,v4?)) = ff(gg,pp)*fi(v1,v2,v3,v4)*fj(p1,p2);
id ff(gg,pp)*fi(?args2) = distrib_(+1,2,gg,pp,?args2);
id ff(gg,gg)*fi(v1?,v2?,v3?,v4?) = 1/2*distrib_(+1,2,gg,gg,v1,v2,v3,v4);
id pp(v1?,v2?)*fj(p1?,p2?) = p1(v1)*p2(v2);
id gg(v1?,v2?) = d_(v1,v2);
.sort

b TI;
.sort
Keep Brackets;
* rank 5 (by Seth, these should be checked by another code that does rank 5)

id TI(n1?,k(v1?),k(v2?),k(v3?),k(v4?),k(v5?), p1?!{k}, ?args) =
    sum_(c1,1,n1-1, 
        sum_(c2,1,n1-1,
            sum_(c3,1,n1-1,
                sum_(c4,1,n1-1,
                    sum_(c5,1,n1-1,
                         f(c1, k(v1), p1, ?args)
                        *f(c2, k(v2), p1, ?args)
                        *f(c3, k(v3), p1, ?args)
                        *f(c4, k(v4), p1, ?args)
                        *f(c5, k(v5), p1, ?args)
                        *TI(n1,sf(c1,c2,c3,c4,c5), p1, ?args)
                    )
                )
            )
        )
    ) + sum_(c1,1,n1-1,
            sum_(c2,1,n1-1,
                sum_(c3,1,n1-1,
                    ff(ddd(v1,v2), f(c1,k(v3),p1,?args), 
                                   f(c2,k(v4),p1,?args),
                                   f(c3,k(v5),p1,?args))
                   *TI(n1,0,0,sf(c1,c2,c3),p1,?args)
                )
            )
        )
    + sum_(c1,1,n1-1,
          ff(ddd(v1,v2), ddd(v3,v4), f(c1,k(v5),p1,?args))
         *TI(n1,0,0,0,0,c1,p1,?args)
      )
;

repeat id once f(n1?{,>1}, k(v1?), p1?!{k}, p2?, ?args) = f(n1-1, k(v1), p1+p2, ?args);
id f(1,k(v1?), p1?!{k}, ?args) = p1(v1);

argument;
    repeat id once f(n1?{,>1}, k(v1?), p1?!{k}, p2?, ?args) = f(n1-1, k(v1), p1 + p2, ?args);
    id f(1,k(v1?), p1?!{k}, ?args) = f(p1,v1);
endargument;

id ff(ddd(v1?,v2?),f(p1?,v3?), f(p2?,v4?), f(p3?,v5?)) = ff(gg,ppp)*fi(v1,v2,v3,v4,v5)*fj(p1,p2,p3);
id ff(ddd(v1?,v2?),ddd(v3?,v4?),f(p1?,v5?)) = ff(gg,gg,pp)*fi(v1,v2,v3,v4,v5)*fj(p1);
id ff(gg,ppp)*fi(?args2) = distrib_(+1,2,gg,ppp,?args2);
id ff(gg,gg,pp)*fi(?args2) = distrib_(+1,4,gggg,pp,?args2);
id ppp(v1?,v2?,v3?)*fj(p1?,p2?,p3?) = p1(v1)*p2(v2)*p3(v3);
id gggg(v1?,v2?,v3?,v4?) = 1/2*distrib_(+1,2,gg,gg,v1,v2,v3,v4);
id pp(v1?)*fj(p1?) = p1(v1);
id gg(v1?,v2?) = d_(v1,v2);
.sort

b TI;
.sort
Keep Brackets;
argument;
    symmetrize sf;
endargument;

b TI;
.sort
Keep Brackets;
id TI(n1?, sf(?args), ?args2) = TI(n1, ?args,?args2);
id TI(n1?, n2?, n3?, sf(?args), ?args2) = TI(n1, n2, n3, ?args,?args2);
.sort

#endprocedure
