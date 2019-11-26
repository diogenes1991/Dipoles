#procedure redefineFields(CONJ,LOOPS)

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
s ccLp;
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
** Change Seth's color correlator insertion code to consider all legs. Currently
** hard-coded for 6 possible external partons. 
*#if `CONJ' != 1
** Insert place holders for color correlations. Currently the sum of all choices 
** that can be made for up to 6 particles. We will not consider diagonal terms, 
** and we use ccLp to keep track of which matches have actually been found in 
** the process at hand. 
*  multiply ccLp*ccLp*
*           ( ccL1*ccL2 +
*             ccL1*ccL3 +
*             ccL1*ccL4 +
*             ccL1*ccL5 +
*             ccL1*ccL6 +
*             ccL2*ccL3 +
*             ccL2*ccL4 +
*             ccL2*ccL5 +
*             ccL2*ccL6 +
*             ccL3*ccL4 +
*             ccL3*ccL5 +
*             ccL3*ccL6 +
*             ccL4*ccL5 +
*             ccL4*ccL6 +
*             ccL5*ccL6 );
*  #do i=1,6
*** Insert appropriate correlators. The adjoint index cOli9e is reserved for the 
*** correlating gluon's adjoint color index.
*** initial quark
*    id QQ(U,p`i',?args,cOli1?)*cc(p`i',cOli1?)*ccLp          *ccL`i' 
*     = QQ(U,p`i',?args,cOli1) *(-cOlT(cOli`i'e,cOli1,cOli9e))*ccL`i';
*** final quark
*    id QQ(Bar(U),p`i',?args,cOli1?)*cc(p`i',cOli1?)*ccLp          *ccL`i' 
*     = QQ(Bar(U),p`i',?args,cOli1) *(+cOlT(cOli1,cOli`i'e,cOli9e))*ccL`i';
*** initial antiquark
*    id QQ(Bar(V),p`i',?args,cOli1?)*cc(p`i',cOli1?)*ccLp          *ccL`i'
*     = QQ(Bar(V),p`i',?args,cOli1) *(+cOlT(cOli1,cOli`i'e,cOli9e))*ccL`i';
*** final antiquark
*    id QQ(V,p`i',?args,cOli1?)*cc(p`i',cOli1?)*ccLp          *ccL`i' 
*     = QQ(V,p`i',?args,cOli1) *(-cOlT(cOli`i'e,cOli1,cOli9e))*ccL`i';
*** gluon
*    id VB(n1?,?args,p`i',cOli1?)*cc(p`i',cOli1?)*ccLp          *ccL`i'
*     = VB(n1,?args,p`i',cOli1)  *(+cOlf(cOli1,cOli9e,cOli`i'e))*ccL`i';        
*** In each term (i.e. color correlated Born ME) only two legs are correlated, 
*** for all other legs in each term we still need to do the usual replacement.
*    id cc(p`i', cOli1?) 
*     = deltaN(cOli1, cOli`i'e);
*  #enddo
*** Remove those correlator terms that are not possible in the process at hand,
*** i.e. those that did not get matched.
*  id ccLp = 0;
*#else
*  #do i=1,6
*    id cc(p`i', cOli1?) 
*     = deltaN(cOli1, cOli`i'e);
*  #enddo
*#endif

** Dio 20191114:
** Color correlated amplitudes. We use ccL0 to track the Born, while ccLi*ccLj 
** track the <TiTj> color correlated Borns.
** CR 20191120:
** Modify the reading direction of the color indices in the quark and antiquark
** color collerators. We always read against the fermion-flow arrows in Feynman 
** diagrams, i.e. for initial quarks / final antiquarks, where the fermion flow 
** arrows point inwards, in the T matrix the inner index comes first, while for 
** initial antiquarks / final quarks, with the fermion flow arrow pointing out-
** wards, the inner index comes second. Directly comment the old ones out below
** in favor of the new ones. For gluons nothing changes for now.  As T matrices 
** are hermitian, I am still confused about the signs in the way Catani-Seymour 
** (CS) define color correlations, but it seems to work (at least when checking 
** against in-house routines and external codes which supposedly also use color
** correlations ala CS). Also only transposing, not hermitian conjugating, i.e. 
** only swapping the two non-adjoint indices, does not result in minus signs in 
** all T matrices.
** Update on the gluon correlator. By trial and error, wrt Recola but also what 
** analytic formulas say (appendix A in CS) we also need to swap the indices in 
** question in the gluon correlator f.
#do i=1,6
  #if `colorCorr' == 1 && `CONJ' != 1 && `LOOPS' == 0
** initial quark
    id QQ(U,p`i',?args,cOli1?)*cc(p`i',cOli1?)
    = QQ(U,p`i',?args,cOli1)*
*      (ccL0*deltaN(cOli1, cOli`i'e) - ccLp*ccL`i'*cOlT(cOli`i'e,cOli1,cOli9e));
      (ccL0*deltaN(cOli1, cOli`i'e) - ccLp*ccL`i'*cOlT(cOli1,cOli`i'e,cOli9e));
** initial antiquark
    id QQ(Bar(V),p`i',?args,cOli1?)*cc(p`i',cOli1?)
    =  QQ(Bar(V),p`i',?args,cOli1)*
*       (ccL0*deltaN(cOli1, cOli`i'e) + ccLp*ccL`i'*cOlT(cOli1,cOli`i'e,cOli9e));
       (ccL0*deltaN(cOli1, cOli`i'e) + ccLp*ccL`i'*cOlT(cOli`i'e,cOli1,cOli9e));
** final quark
   id QQ(Bar(U),p`i',?args,cOli1?)*cc(p`i',cOli1?)
   =  QQ(Bar(U),p`i',?args,cOli1)*
*      (ccL0*deltaN(cOli1, cOli`i'e) + ccLp*ccL`i'*cOlT(cOli1,cOli`i'e,cOli9e));
      (ccL0*deltaN(cOli1, cOli`i'e) + ccLp*ccL`i'*cOlT(cOli`i'e,cOli1,cOli9e));
** final antiquark
    id QQ(V,p`i',?args,cOli1?)*cc(p`i',cOli1?)
    = QQ(V,p`i',?args,cOli1)*
*      (ccL0*deltaN(cOli1, cOli`i'e) - ccLp*ccL`i'*cOlT(cOli`i'e,cOli1,cOli9e));
      (ccL0*deltaN(cOli1, cOli`i'e) - ccLp*ccL`i'*cOlT(cOli1,cOli`i'e,cOli9e));
** gluon
    id VB(n1?,?args,p`i',cOli1?)*cc(p`i',cOli1?)
    = VB(n1,?args,p`i',cOli1)*
*      (ccL0*deltaN(cOli1, cOli`i'e) + ccLp*ccL`i'*i_*cOlf(cOli1,cOli9e,cOli`i'e));
      (ccL0*deltaN(cOli1, cOli`i'e) + ccLp*ccL`i'*i_*cOlf(cOli`i'e,cOli9e,cOli1));
** we have to throw away everything with more powers than two of the correlators
    if(count(ccLp,1)>2)discard;
** to keep track of the Born and of the non-correlated legs, we don't need 
** multiple powers of ccL0
    if(count(ccL0,1)==2) id ccL0^2 = ccL0;
  #else
    id cc(p`i', cOli1?) = deltaN(cOli1, cOli`i'e);
  #endif
#enddo
** Disentangle the Born from the color-correlated Borns. If we kept ccL0 also in 
** those terms, then we would not only return the Born when asking for ccL0 = 1.
#do i=1,6
  id ccL`i'*ccL0 = ccL`i';
#enddo 
** Discard terms with only one color correlator and throw away power counter.
if(count(ccLp,1)==1) discard;
id ccLp = 1;

b deltaN, cOlT, cOlf, delta, cOlOne, VB;
.sort 
Keep Brackets;

** Comment by Chris on 20161107:
** With the Keep Bracket statement above, the following statements are supposed 
** to be only applied to the contents outside the bracket. Remember that b X,Y; 
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
