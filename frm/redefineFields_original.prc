#procedure redefineFields(CONJ)

cf cc;

* TODO: replace below by auto-generated from model file

id Uu(?args) = QQ(U, ?args);
id Vu(?args) = QQ(V, ?args);
id Ud(?args) = QQ(U, ?args);
id Vd(?args) = QQ(V, ?args);
id Ub(?args) = QQ(U, ?args);
id Vb(?args) = QQ(V, ?args);
id Ut(?args) = QQ(U, ?args);
id Vt(?args) = QQ(V, ?args);
id Ul(?args) = QQ(U, ?args);
id Vl(?args) = QQ(V, ?args);
id Uq(?args) = QQ(U, ?args);
id Vq(?args) = QQ(V, ?args);
id UuBar(?args) = QQ(Bar(U), ?args);
id VuBar(?args) = QQ(Bar(V), ?args);
id UdBar(?args) = QQ(Bar(U), ?args);
id VdBar(?args) = QQ(Bar(V), ?args);
id UbBar(?args) = QQ(Bar(U), ?args);
id VbBar(?args) = QQ(Bar(V), ?args);
id UtBar(?args) = QQ(Bar(U), ?args);
id VtBar(?args) = QQ(Bar(V), ?args);
id UlBar(?args) = QQ(Bar(U), ?args);
id VlBar(?args) = QQ(Bar(V), ?args);
id UqBar(?args) = QQ(Bar(U), ?args);
id VqBar(?args) = QQ(Bar(V), ?args);

id EpsG(?args) = VB(fEpsG, ?args);
id EpsGStar(?args) = VB(fEpsGStar, ?args);
id EpsA(?args) = VB(fEpsA, ?args);
id EpsAStar(?args) = VB(fEpsAStar, ?args);
id EpsZ(?args) = VB(fEpsZ, ?args);
id EpsZStar(?args) = VB(fEpsZStar, ?args);
id EpsW(?args) = VB(fEpsW, ?args);
id EpsWStar(?args) = VB(fEpsWStar, ?args);

id QQ(n1?, p1?, ?args, cOli1?) = QQ(n1, p1, ?args, cOli1)*cc(p1,cOli1);
id VB(n1?, ?args, p1?, cOli1?) = VB(n1, ?args, p1, cOli1)*cc(p1,cOli1);

** Comment by Chris on 20161107:
** The above id for VB(n1?, ?args, p1?, cOlil1?) is rather misleading. The VB in
** the above refers obviously to a gluon. However, gluons don't have fundamental
** color indices with dimension N, which one could be led to believe from the id
** cc(...) = deltaN(...) below, but adjoint color indices with dim N^2-1. Though
** it seems that deltaN is used only locally, and that, when used in conjunction
** with possible adjoint indices (see next repeat statement) it is done correct.

** Comment by CR on 20191106 (addition to comment above):
** What we have above is, that color indices also of a VB (with momentum pX and 
** color index cOliY) gets a cc(pX,cOliY) assigned. Further belowi cc(pX,cOliY) 
** is then identified with deltaN(cOliY,cOliXe). Using deltaN for this is a bit 
** misleading, as the N in deltaN suggests this to be a delta for color indices 
** in the fundamental representation only (and deltaN may actually have certain 
** properties that rely on that fact). However in here deltaN is simply used as 
** a mere delta to identify two color indices, irregardless whether they are in
** the fundamental or adjoint representation, in the repeat block at the bottom
** and then set to one. No other contractions of indices, which would lead to a
** definite number, are performed with deltaN in here.

** Comment by CR on 20191106:
** cOlT(i,j,a) takes two fundamental indices, then one or more adjoint indices.
** Have a look, for example, at the identities in the files NLOXSimplifyColor.h 
** or NLOXSimplifyColorMatrices.prc, etc. 
** Unfortunately, we do not always pay particular attention to denoting the in-
** dices accordingly. We should, however, use some fundamental indices with di-
** mension N for the first two indices and adjoint indices with dimension N^2-1 
** for the other (say some cOli and cOla; we could even distinguish fundamental 
** from anti-fundamental, say cOli and cOlj).

b cc,QQ,VB;
.sort
cf deltaN;
s ccLp,ccRp;
dimension cOlNA;
Keep Brackets;
*#do i=1,6
** Insert color correlator matrix if appropriate.
** cOli9e is reserved for the extra particle's color index;
** change when hard-coded bounds are removed
*    #if `CONJ' != 1 && ( ( `i' == `insertleg1' ) || ( `i' == `insertleg2' ) )
** initial quark
*        id QQ(U,p`i',?args,cOli1?)*cc(p`i',cOli1?) = -QQ(U,p`i',?args,cOli1)*cOlT(cOli`i'e,cOli1,cOli9e);
** final quark
*        id QQ(Bar(U),p`i',?args,cOli1?)*cc(p`i',cOli1?) = QQ(Bar(U),p`i',?args,cOli1)*cOlT(cOli1,cOli`i'e,cOli9e);
** initial antiquark
*        id QQ(Bar(V),p`i',?args,cOli1?)*cc(p`i',cOli1?) = QQ(Bar(V),p`i',?args,cOli1)*cOlT(cOli1,cOli`i'e,cOli9e);
** final antiquark
*        id QQ(V,p`i',?args,cOli1?)*cc(p`i',cOli1?) = -QQ(V,p`i',?args,cOli1)*cOlT(cOli`i'e,cOli1,cOli9e);
** gluon
*        id VB(n1?,?args,p`i',cOli1?)*cc(p`i',cOli1?) = VB(n1,?args,p`i',cOli1)*cOlf(cOli1,cOli9e,cOli`i'e);        
*    #else
*        id cc(p`i', cOli1?) = deltaN(cOli1, cOli`i'e);
*    #endif
*#enddo
** CR 20191108:
** Change Seth's color correlator insertion code to consider all legs. Use two 
** different identifiers ccR and ccL for now, to separate the insertions right 
** (conj) and left (non-conj) -- ccR could also be non-conj and ccL conj. Also
** introduced two new symbols ccLp and ccRp above, locally, to count powers of 
** insertions.
#do i=1,6
** Insert color correlator matrix if appropriate.
** cOli9e is reserved for the extra particle's color index;
** change when hard-coded bounds are removed
  #if `CONJ' != 1
** initial quark
    id QQ(U,p`i',?args,cOli1?)*
       cc(p`i',cOli1?)
     = QQ(U,p`i',?args,cOli1)*
       cc(p`i',cOli1)*
       (1-cOlT(cOli`i'e,cOli1,cOli9e)*ccL`i'*ccLp)*
       (1-cOlT(cOli`i'e,cOli1,cOli9e)*ccR`i'*ccRp);
** final quark
    id QQ(Bar(U),p`i',?args,cOli1?)*
       cc(p`i',cOli1?)
     = QQ(Bar(U),p`i',?args,cOli1)*
       cc(p`i',cOli1)*
       (1+cOlT(cOli1,cOli`i'e,cOli9e)*ccL`i'*ccLp)*
       (1+cOlT(cOli1,cOli`i'e,cOli9e)*ccR`i'*ccRp);
** initial antiquark
    id QQ(Bar(V),p`i',?args,cOli1?)*
       cc(p`i',cOli1?)
     = QQ(Bar(V),p`i',?args,cOli1)*
       cc(p`i',cOli1)*
       (1+cOlT(cOli1,cOli`i'e,cOli9e)*ccL`i'*ccLp)*
       (1+cOlT(cOli1,cOli`i'e,cOli9e)*ccR`i'*ccRp);
** final antiquark
    id QQ(V,p`i',?args,cOli1?)*
       cc(p`i',cOli1?)
     = QQ(V,p`i',?args,cOli1)*
       cc(p`i',cOli1)*
       (1-cOlT(cOli`i'e,cOli1,cOli9e)*ccL`i'*ccLp)*
       (1-cOlT(cOli`i'e,cOli1,cOli9e)*ccR`i'*ccRp);
** gluon
    id VB(n1?,?args,p`i',cOli1?)*
       cc(p`i',cOli1?)
     = VB(n1,?args,p`i',cOli1)*
       cc(p`i',cOli1)*
       (1+cOlf(cOli1,cOli9e,cOli`i'e)*ccL`i'*ccLp)*
       (1+cOlf(cOli1,cOli9e,cOli`i'e)*ccR`i'*ccRp);
  #else
    id cc(p`i', cOli1?) 
     = deltaN(cOli1, cOli`i'e);
  #endif
#enddo
** Only keep those terms that have one ccL AND one ccR. As the above code does 
** not make sense for the diagonal terms of the color correlator matrix, since
** other color indices are also used twice in that case, we remove those also.
#if `CONJ' != 1
** Get rid of those terms that have the wrong number of insertions.
  if(count(ccLp,1)>1) discard;
  if(count(ccRp,1)>1) discard;
  if(match(ccLp)=0) discard;
  if(match(ccRp)=0) discard;
** Set reamining factors that are not needed anymore to 1.
  id ccLp=1;
  id ccRp=1;
  id cc(?args)=1;
** Get rid of the diagonal terms of the color correlator matrix. Actually, also
** get rid of one of the non-diagonal triangles (currently such that only terms 
** with ccL`i'*ccR`k' for i>k remain). 
  #do i=1,6
  id ccL`i'*ccR`i' = 0;
  #do k={`i'+1},6
    id ccL`i'*ccR`k' = 0;
  #enddo 
  #enddo 
#endif
b deltaN, cOlT, cOlf, delta, cOlOne, VB;
.sort 
Keep Brackets;

** Comment by Chris on 20161107:
** With the Keep Bracket statement above, the following statements are supposed 
** to be only applied to the contens outside the bracket. Remember that b X, Y; 
** really means to bracket-out X, Y.

repeat;
  id once cOlT(cOli1?, cOli2?, cOli3?)*deltaN(cOli1?, cOli4?) = cOlT(cOli4, cOli2, cOli3)*deltaN(cOli1, cOli4);
  id once cOlT(cOli2?, cOli1?, cOli3?)*deltaN(cOli1?, cOli4?) = cOlT(cOli2, cOli4, cOli3)*deltaN(cOli1, cOli4);
  id once cOlT(cOli2?, cOli3?, cOli1?)*deltaN(cOli1?, cOli4?) = cOlT(cOli2, cOli3, cOli4)*deltaN(cOli1, cOli4);
  id once cOlf(cOli1?, cOli2?, cOli3?)*deltaN(cOli1?, cOli4?) = cOlf(cOli4, cOli2, cOli3)*deltaN(cOli1, cOli4);
  id once cOlf(cOli2?, cOli1?, cOli3?)*deltaN(cOli1?, cOli4?) = cOlf(cOli2, cOli4, cOli3)*deltaN(cOli1, cOli4);
  id once cOlf(cOli2?, cOli3?, cOli1?)*deltaN(cOli1?, cOli4?) = cOlf(cOli2, cOli3, cOli4)*deltaN(cOli1, cOli4);
  id once deltaN(cOli1?, cOli2?)*delta(cOli1?, cOli3?) = deltaN(cOli1, cOli2)*delta(cOli2, cOli3);
  id once deltaN(cOli1?, cOli2?)*VB(?args, cOli1?) = deltaN(cOli1, cOli2)*VB(?args, cOli2);
endrepeat;
b deltaN;
.sort
id deltaN(cOli1?, cOli2?) = 1;

b cOlT, cOlf;
.sort

b VB;
.sort
dimension d;

#endprocedure
