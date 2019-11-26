#procedure NLOXSimplifyColorMatrices

** Meant as alternative to colorMatrices.prc.

** Note: In the following we rather rely on the pre-defined functions cOlf, cOlT
** and cOlTr rather than defining them ourselves, since this routine is meant to
** be used in evaluateDiags,  at which stage these objects have been defined via
** #include color.h in declarations.h by
** color.h:Tensor cOlfp,cOlff(cyclic),cOlT,cOlTr(cyclic),cOlTt(cyclic),cOlf3;
** color.h:Tensor cOlf(antisymmetric),cOlf1(antisymmetric);

.sort
** Just so I don't have to rewrite everything. However here cOla is not assigned 
** a definite dimension of NA.
auto i cOla;
** It seems I cannot directly use cOlTr, being defined as cyclic tensor, in con-
** junction with argument field wildcards. However, if I intermittently use ano-
** ther function and then identify this one with cOlTr, it seems to work.
cf cOlTra;

#message NLOXSimplifyColorMatrices
b cc;
.sort:cc bracket;
id cc(?args) = cc(cc(?args));

b cOlT, cOlf, delta, cOlOne;
.sort: cc wrapping;

** On the diagram level, except for a closed quark loop with no gluons attached,
** there should be no full contraction of color indices in Kronecker deltas that
** would lead to CA=NC or NA=NC^2-1,  and hence the following should generically 
** be applicable. In the following we explicitely use a normalization of 1/2 and
** NC=3.

id delta(cOli1?, cOli2?) = d_(cOli1, cOli2);
.sort

** The following numbers in parantheses refer to the equations listed in the be-
** ginning of NLOXSimplifyColor.h.

** Operations on color matrices in the adjoint representation first
** (1)
id cOlf(cOli1?,cOli2?,cOli3?) = -i_*2*(cOlTra(cOli1,cOli2,cOli3)-cOlTra(cOli1,cOli3,cOli2));
.sort
** Operations on traces that originate from the adjoint color matrices.
repeat;
** (3)
  id cOlTra(?X,cOli1?,?Y)*cOlTra(?A,cOli1?,?B)=1/2*(cOlTra(?Y,?X,?B,?A)-1/3*cOlTra(?Y,?X)*cOlTra(?B,?A));
** (4)
  id cOlTra(?A,cOli1?,?X,cOli1?,?Y)=1/2*(cOlTra(?X)*cOlTra(?Y,?A)-1/3*cOlTra(?X,?Y,?A));
** (5) and (6)
*  id cOlTra(cOli1?,cOli2?)=1/2*d_(cOli1,cOli2);
  id cOlTra(cOli1?,cOli2?)=1/2*delta(cOli1,cOli2);
  id cOlTra(cOli1?)=0;
  id cOlTra()=3;
** One more needed here
*  id cOlT(cOli1?,cOli2?,cOli3?)*cOlTra(?X,cOli3?,?Y) = 1/2*(cOlT(cOli1,cOli2,?Y,?X)-1/3*cOlTra(?Y,?X)*d_(cOli1,cOli2));
*  id cOlT(cOli1?,cOli2?,?A1,cOli5?,?A2)*cOlT(cOli3?,cOli4?,?B1,cOli5?,?B2) = 1/2*(cOlT(cOli1,cOli4,?A1,?B2)*cOlT(cOli3,cOli2,?B1,?A2)-1/3*cOlT(cOli1,cOli2,?A1,?A2)*cOlT(cOli3,cOli4,?B1,?B2));
*  id cOlT(cOli1?,cOli2?) = d_(cOli1,cOli2);
  id cOlT(cOli1?,cOli2?,cOli3?)*cOlTra(?X,cOli3?,?Y) = 1/2*(cOlT(cOli1,cOli2,?Y,?X)-1/3*cOlTra(?Y,?X)*delta3(cOli1,cOli2));
  id cOlT(cOli1?,cOli2?,?A1,cOli5?,?A2)*cOlT(cOli3?,cOli4?,?B1,cOli5?,?B2) = 1/2*(cOlT(cOli1,cOli4,?A1,?B2)*cOlT(cOli3,cOli2,?B1,?A2)-1/3*cOlT(cOli1,cOli2,?A1,?A2)*cOlT(cOli3,cOli4,?B1,?B2));
  id cOlT(cOli1?,cOli2?) = delta3(cOli1,cOli2);
endrepeat;
.sort

** Operations on objects with explicit fundamental color indices
repeat;
** (2)
*  id cOlT(cOli1?,cOli2?,cOla1?)*cOlT(cOli3?,cOli4?,cOla1?) = 1/2*(d_(cOli1,cOli4)*d_(cOli3,cOli2)-1/3*d_(cOli1,cOli2)*d_(cOli3,cOli4));
*  id cOlT(cOli1?,cOli2?,cOli5?)*cOlT(cOli3?,cOli4?,cOli5?) = 1/2*(d_(cOli1,cOli4)*d_(cOli3,cOli2)-1/3*d_(cOli1,cOli2)*d_(cOli3,cOli4));
  id cOlT(cOli1?,cOli2?,cOla1?)*cOlT(cOli3?,cOli4?,cOla1?) = 1/2*(delta3(cOli1,cOli4)*delta3(cOli3,cOli2)-1/3*delta3(cOli1,cOli2)*delta3(cOli3,cOli4));
  id cOlT(cOli1?,cOli2?,cOli5?)*cOlT(cOli3?,cOli4?,cOli5?) = 1/2*(delta3(cOli1,cOli4)*delta3(cOli3,cOli2)-1/3*delta3(cOli1,cOli2)*delta3(cOli3,cOli4));
** (5)
*  id cOlT(cOli1?,cOli2?,cOla1?)*cOlT(cOli2?,cOli1?,cOla2?) = 1/2*d_(cOla1,cOla2);
*  id cOlT(cOli1?,cOli2?,cOli3?)*cOlT(cOli2?,cOli1?,cOli4?) = 1/2*d_(cOli3,cOli4);
  id cOlT(cOli1?,cOli2?,cOla1?)*cOlT(cOli2?,cOli1?,cOla2?) = 1/2*delta(cOla1,cOla2);
  id cOlT(cOli1?,cOli2?,cOli3?)*cOlT(cOli2?,cOli1?,cOli4?) = 1/2*delta(cOli3,cOli4);
  id cOlT(cOli1?,cOli1?,cOla1?) = 0;
  id cOlT(cOli1?,cOli1?,cOli2?) = 0;
**
  id delta(cOli1?,cOli2) = d_(cOli1,cOli2);
endrepeat;
.sort

** Identify left over traces
repeat;
  id cOlT(cOli1?,cOli2?,cOla1?)*cOlT(cOli2?,cOli1?,cOla2?) = cOlTra(cOla1,cOla2);
  id cOlT(cOli1?,cOli2?,cOla1?)*cOlT(cOli2?,cOli3?,cOla2?)*cOlT(cOli3?,cOli1?,cOla3?) = cOlTra(cOla1,cOla2,cOla3);
  id cOlT(cOli1?,cOli2?,cOla1?)*cOlT(cOli2?,cOli3?,cOla2?)*cOlT(cOli3?,cOli4?,cOla3?)*cOlT(cOli4?,cOli1?,cOla4?) = cOlTra(cOla1,cOla2,cOla3,cOla4);
  id cOlT(cOli1?,cOli2?,cOla1?)*cOlT(cOli2?,cOli3?,cOla2?)*cOlT(cOli3?,cOli4?,cOla3?)*cOlT(cOli4?,cOli5?,cOla4?)*cOlT(cOli5?,cOli1?,cOla5?) = cOlTra(cOla1,cOla2,cOla3,cOla4,cOla5);
  #do i = 5, 30
    id cOlT(cOli1?,cOli2?,cOla1?)*<cOlT(cOli2?,cOli3?,cOla2?)>*...*<cOlT(cOli`i'?,cOli{`i'+1}?,cOla`i'?)>*cOlT(cOli{`i'+1}?,cOli1?,cOla{`i'+1}?) = cOlTra(cOla1,cOla2,...,cOla`i',cOla{`i'+1});
  #enddo
  id cOlT(cOli1?,cOli2?,cOli3?)*cOlT(cOli2?,cOli1?,cOli4?) = cOlTra(cOli3,cOli4);
  id cOlT(cOli1?,cOli2?,cOli4?)*cOlT(cOli2?,cOli3?,cOli5?)*cOlT(cOli3?,cOli1?,cOli6?) = cOlTra(cOli4,cOli5,cOli6);
  id cOlT(cOli1?,cOli2?,cOli5?)*cOlT(cOli2?,cOli3?,cOli6?)*cOlT(cOli3?,cOli4?,cOli7?)*cOlT(cOli4?,cOli1?,cOli8?) = cOlTra(cOli5,cOli6,cOli7,cOli8);
  id cOlT(cOli1?,cOli2?,cOli6?)*cOlT(cOli2?,cOli3?,cOli7?)*cOlT(cOli3?,cOli4?,cOli8?)*cOlT(cOli4?,cOli5?,cOli9?)*cOlT(cOli5?,cOli1?,cOli10?) = cOlTra(cOli6,cOli7,cOli8,cOli9,cOli10);
  #do i = 5, 30
    id cOlT(cOli1?,cOli2?,cOli{`i'+2}?)*<cOlT(cOli2?,cOli3?,cOli{`i'+3}?)>*...*<cOlT(cOli`i'?,cOli{`i'+1}?,cOli{`i'+`i'+1}?)>*cOlT(cOli{`i'+1}?,cOli1?,cOli{`i'+1+`i'+1}?) = cOlTra(cOli{`i'+2},cOli{`i'+3},...,cOli{`i'+`i'+1},cOli{`i'+1+`i'+1});
  #enddo
endrepeat;
.sort

** Operations on left over traces
repeat;
** (3)
  id cOlTra(?X,cOli1?,?Y)*cOlTra(?A,cOli1?,?B)=1/2*(cOlTra(?Y,?X,?B,?A)-1/3*cOlTra(?Y,?X)*cOlTra(?B,?A));
** (4)
  id cOlTra(?A,cOli1?,?X,cOli1?,?Y)=1/2*(cOlTra(?X)*cOlTra(?Y,?A)-1/3*cOlTra(?X,?Y,?A));
** (5) and (6)
*  id cOlTra(cOli1?,cOli2?)=1/2*d_(cOli1,cOli2);
  id cOlTra(cOli1?,cOli2?)=1/2*delta(cOli1,cOli2);
  id cOlTra(cOli1?)=0;
  id cOlTra()=3;
**
  id delta(cOli1?,cOli2) = d_(cOli1,cOli2);
endrepeat;
.sort

id cOlTra(?X)=cOlTr(?X);
*id d_(cOli1?,cOli2?) = delta(cOli1,cOli2);
.sort
id delta3(cOli1?,cOli2?) = delta(cOli1,cOli2);
.sort
id delta(cOli1?,cOli2?) = d_(cOli1,cOli2);
.sort
id d_(cOli1?,cOli2?) = delta(cOli1,cOli2);
.sort
id cOlT(cOli1?,cOli2?,cOli5?)*cOlT(cOli3?,cOli4?,cOli5?) = 1/2*(delta3(cOli1,cOli4)*delta3(cOli3,cOli2)-1/3*delta3(cOli1,cOli2)*delta3(cOli3,cOli4));
.sort
id delta3(cOli1?,cOli2?) = delta(cOli1,cOli2);
.sort
id delta(cOli1?,cOli2?) = d_(cOli1,cOli2);
.sort
id d_(cOli1?,cOli2?) = delta(cOli1,cOli2);
.sort

print +s;
.sort

b cc;
.sort:final cc bracket;
dimension d;
id cc(cc(?args)) = cc(?args);

#endprocedure

