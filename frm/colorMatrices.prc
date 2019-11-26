#procedure colorMatrices
#message colormatrices
b cc;
.sort:cc bracket;
dimension cOlNA;
id cc(?args) = cc(cc(?args));

b cOlT, cOlf, delta, cOlOne;
.sort:cc wrapping;
s ttt;
Keep Brackets;

repeat id delta(cOli1?, cOli2?)*delta(cOli2?, cOli3?) = delta(cOli1, cOli3);
id delta(cOli1?, cOli1?) = 3;

id delta(cOli1?, cOli2?) = d_(cOli1, cOli2);

*
* WARNING: delta3 dependent stuff produces wrong stuff
*
*

#do l=1,4
    b cOlT, cOlf, delta, cOlOne, delta3;
    .sort:color loop;
    s ttt;
    Keep Brackets;
*    id cOlf(cOli3?, cOli2?, cOli1?)*cOlT(cOli4?, cOli5?, cOli3?) =
*         -i_*(cOlT(cOli4, cOli{`l'+5}, cOli1)*cOlT(cOli{`l'+5}, cOli5, cOli2) - cOlT(cOli4, cOli{`l'+5}, cOli2)*cOlT(cOli{`l'+5}, cOli5, cOli1));
    repeat;
       id cOlT(cOli1?, cOli2?, cOli3?)*cOlT(cOli2?, cOli1?, cOli4?) = 1/2*d_(cOli3, cOli4);
       id cOlT(cOli1?, cOli1?, cOli2?) = 0;
       id cOlf(cOli1?, cOli2?, cOli3?)*cOlf(cOli4?, cOli2?, cOli3?) = CA*d_(cOli1, cOli4);
       id cOlf(cOli1?, cOli2?, cOli3?)*cOlT(cOli4?, cOli5?, cOli2?)*cOlT(cOli5?, cOli6?, cOli3?) = 
            1/2*i_*CA*cOlT(cOli4, cOli6, cOli1);
*       id cOlT(cOli1?, cOli2?, cOli3?)*cOlT(cOli4?, cOli5?, cOli3?) = 
*           1/2*(delta3(cOli1, cOli5)*delta3(cOli2, cOli4) - 1/3*delta3(cOli1, cOli2)*delta3(cOli4, cOli5));
       id delta3(cOli1?, cOli2?)*delta3(cOli2?,cOli3?) = delta3(cOli1, cOli3);
       id delta3(cOli1?, cOli1?) = 3;
       id delta3(cOli2?, cOli1?)*cOlT(cOli1?, cOli3?, cOli4?) = cOlT(cOli2, cOli3, cOli4);
       id delta3(cOli2?, cOli1?)*cOlT(cOli3?, cOli1?, cOli4?) = cOlT(cOli3, cOli2, cOli4);
    endrepeat;
#enddo

b d_, cOlOne, delta3, cOlT;
.sort:color d_ bracket;
Keep Brackets;
id d_(cOli1?, cOli2?) = delta(cOli1, cOli2);

repeat;
    id delta3(cOli1?, cOli2?)*delta3(cOli2?,cOli3?) = delta3(cOli1, cOli3);
    id delta3(cOli1?, cOli1?) = 3;
    id delta3(cOli2?, cOli1?)*cOlT(cOli1?, cOli3?, cOli4?) = cOlT(cOli2, cOli3, cOli4);
    id delta3(cOli2?, cOli1?)*cOlT(cOli3?, cOli1?, cOli4?) = cOlT(cOli3, cOli2, cOli4);
endrepeat;


b cc;
.sort:final cc bracket;
dimension d;
id cc(cc(?args)) = cc(?args);

#endprocedure
