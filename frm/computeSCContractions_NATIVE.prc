#procedure computeSCContractions

#define SPINLINES "3"
* Note: NUMERICALTRACES should be left to off for now, because it is not fully
* implemented. Return to development once we consider doing numerical
* polarizations as well.
#define NUMERICALTRACES "0"
#include config.inc
#include `leftDiagDefFile' # global
#call eliminateMomentum
.sort

**** Notes on implementation of certain schemes
* "IR" schemes refer to the strictly original HV scheme where all loop boson
* helicities are d-dimensional.

* Non-"IR" schemes have had bosons connecting different fermion lines
* in the amplitude truncated to 4 dimensions--this matches other codes, but
* drops O(eps) contributions of a purely IR nature (hence the temporary name).

* In the code below, we assume that for HV schemes of any type all contractions
* of d-dimensional indices possible before interferences have been done, so be
* careful if modifiying other scripts in a way such that it would prevent
* a repeated index on the same fermion line from being contracted. This
* assumption was made for speed and simplicity.

* "HV" is the presently working but somewhat messy implementation of HV,
* especially its IR version.

#if (`leftNLoops' == 0) || (`dimScheme' == HV)
    dimension 4;
#else * make sure we are still in d
    dimension d;
#endif

id dummySME = 1;

#if `leftNLoops' == 1
b SCC;
.sort

#if (`dimScheme' == HVIR)
* Isolate momenta so they are not matched
repeat id SCC(?args1,G(n0l?,p1?),?args2) = SCC(?args1,G(n0l,acc(p1)),?args2);
* Protect 4 dim indices
repeat id SCC(?args1,G(n0l?,v1?ve),?args2) = SCC(?args1,G(n0l,acc(v1)),?args2);
*#endif

*#if `dimScheme' == HVIR
* 4-dim indices running between fermion lines have been tagged;
* all others between fermion lines should be labeled as a d-dimensional index
    #do i=0,9
        if( match(SCC(n0l?nl,?args1,G(n0l?,v0l?!vdhv$temp),?args2)*
                  SCC(n1l?nl,?args3,G(n1l?,v0l?),?args4)) );
            multiply replace_($temp,vdhv`i');
        elseif( match(SCC(n0l?nl,?args1,G(n0l?,v0l?!vdhv$temp),?args2,G(n0l?,v0l?),?args3)) );
            multiply replace_($temp,vdhv`i');
        endif;
    #enddo

*#endif * `dimScheme' == HVIR

*#if (`dimScheme' == HVIR)
* Return protected indices to normal
repeat id SCC(?args1,G(n0l?,acc(v1?)),?args2) = SCC(?args1,G(n0l,v1),?args2);
repeat id SCC(?args1,G(n0l?,acc(p1?)),?args2) = SCC(?args1,G(n0l,p1),?args2);

* Indices not replaced above are 4-dim; replace with unsummable indices, so
* d-dim metric tensors do not "convert" them to d
    #do i = 0,19
        multiply replace_(v`i'l,vhv`i'l);
        multiply replace_(v`i'r,vhv`i'r);
    #enddo
.sort
#endif * `dimScheme' == HVIR

#endif * `leftNLoops' == 1

#if (`leftNLoops' == 0) || (`dimScheme' == HV)
* All remaining indices should be 4-dim in HV scheme, d-dim in CDR; 
* replace with default dimension
    #do i = 0,19
        sum v`i'l,v`i'r,v`i'el;
    #enddo
.sort
#endif

#call preprocChains

b SCC;
.sort
cf ddf;
Keep Brackets;

#message spin sums
repeat;

* Should be able to handle generic quark spin sum by using mass()

id SCC(dd(n1?,n2?),?args1, QQ(U, p1?, ?args2), QQ(Bar(U), p1?,?args3), ?args4) =
    SCC(dd(n1,n2),?args1, G(n1,p1), ?args4) + mass(p1)*SCC(dd(n1,n2),?args1, G(n1), ?args4);

id SCC(dd(n1?,n2?),?args1, QQ(V, p1?, ?args2), QQ(Bar(V), p1?,?args3), ?args4) =
    SCC(dd(n1,n2),?args1, G(n1,p1), ?args4) - mass(p1)*SCC(dd(n1,n2),?args1, G(n1), ?args4);

endrepeat;

repeat id SCC(dd(n1?,n2?), dd(n3?,n4?), ?args) = SCC(ddf(n1,n2), ddf(n3,n4), ?args);
repeat id SCC(dd(n1?,n2?), ?args) = dd(n1,n2)*SCC(?args);
repeat id SCC(ddf(n1?,n2?), ?args) = ddf(n1,n2)*SCC(?args);

id SCC(?args) = SC(?args);

b SC;
.sort
#message replace chain indices
Keep Brackets;
repeat id once SC(n1?,?args) = n1*SC(?args);
id SC = 1;

b dd, ddf;
.sort
Keep Brackets;
#do i = 0, 10
    id dd(n`i'r, n1?) = dd(n1, n`i'r);
    id ddf(n`i'r, n1?) = ddf(n1, n`i'r);
#enddo
id ddf(n1?, n2?)*ddf(n3?!{,n1}, n4?) = dd(n1, n2)*dd(n1, n4)*dd(n1, n3);

b G, G4, G5, dd;
.sort
f Gc, G5c;
Keep Brackets;
id G4(?args) = G(?args);
repeat;
    id dd(n1?,n2?)*G(n2?,v1?) = Gc(n1,v1)*dd(n1,n2);
    id dd(n1?,n2?)*G(n2?) = Gc(n1)*dd(n1,n2);
    id dd(n1?,n2?)*G5(n2?) = G5c(n1)*dd(n1,n2);
endrepeat;
id dd(n1?,n2?) = 1;

b Gc, G5c;
.sort
Keep Brackets;

id G5c(n1?) = G5(n1);
id Gc(n1?) = G(n1);
id Gc(n1?, v1?) = G(n1, v1);

#call postprocChains

b VB;
.sort

#message >> polarization sums

*id VB(fEpsAStar, v1?, p1?)*VB(fEpsA, v2?, p1?) = -d_(v1,v2);
id VB(fEpsWStar, v1?, p1?)*VB(fEpsW, v2?, p1?) = -d_(v1,v2) + p1(v1)*p1(v2)/(mW*mW);
id VB(fEpsZStar, v1?, p1?)*VB(fEpsZ, v2?, p1?) = -d_(v1,v2) + p1(v1)*p1(v2)/(mZ*mZ);

* Pair of initial state gluons first (unambiguous)
id  VB(fEpsGStar, v1?, p1?{p1,p2})*VB(fEpsG, v2?, p1?)
   *VB(fEpsGStar, v3?, p2?{p1,p2})*VB(fEpsG, v4?, p2?)
 =  (-d_(v1,v2) + 2*DS(p1+p2,0,1)*(p1(v1)*p2(v2) + p1(v2)*p2(v1)))
   *(-d_(v3,v4) + 2*DS(p1+p2,0,1)*(p1(v3)*p2(v4) + p1(v4)*p2(v3)));
id  VB(fEpsAStar, v1?, p1?{p1,p2})*VB(fEpsA, v2?, p1?)
   *VB(fEpsAStar, v3?, p2?{p1,p2})*VB(fEpsA, v4?, p2?)
 =  (-d_(v1,v2) + 2*DS(p1+p2,0,1)*(p1(v1)*p2(v2) + p1(v2)*p2(v1)))
   *(-d_(v3,v4) + 2*DS(p1+p2,0,1)*(p1(v3)*p2(v4) + p1(v4)*p2(v3)));

* Pair of final state gluons (ambiguous if more than two final state gluons)
*
* We use a do loop here in order to avoid the ambiguity that comes with multiple
* possible pairings of gluons. For 3 final-state gluons (p3, p4, p5), FORM could
* arbitrarily pair p3 and p4 for some terms and p4 and p5 for others. This
* breaks gauge invariance. With the do loop, p3 and p4 will be paired for every
* term.
#do i = 3, 5
#do j = 4, 6
id  VB(fEpsGStar, v1?, p`i')*VB(fEpsG, v2?, p`i')
   *VB(fEpsGStar, v3?, p`j')*VB(fEpsG, v4?, p`j')
 =  (-d_(v1,v2) + 2*DS(p`i'+p`j',0,1)*(p`i'(v1)*p`j'(v2) + p`i'(v2)*p`j'(v1)))
   *(-d_(v3,v4) + 2*DS(p`i'+p`j',0,1)*(p`i'(v3)*p`j'(v4) + p`i'(v4)*p`j'(v3)));
id  VB(fEpsAStar, v1?, p`i')*VB(fEpsA, v2?, p`i')
   *VB(fEpsAStar, v3?, p`j')*VB(fEpsA, v4?, p`j')
 =  (-d_(v1,v2) + 2*DS(p`i'+p`j',0,1)*(p`i'(v1)*p`j'(v2) + p`i'(v2)*p`j'(v1)))
   *(-d_(v3,v4) + 2*DS(p`i'+p`j',0,1)*(p`i'(v3)*p`j'(v4) + p`i'(v4)*p`j'(v3)));
#enddo
#enddo

* One initial state gluon, one final state gluon
*
* After identifying intial-state and final-state pairs, there will be at most
* one gluon in the initial state and one in the final state, so no do loop is
* necessary.
*
* Note the different sign for the DS and p1-p2 versus the other two cases.
* Although algebraically equivalent, it seems to be necessary because some DS[]
* with a sum of momenta will not be handled correctly and will end up in the C
* code. The sign difference here prevents that error.
id  VB(fEpsGStar, v1?, p1?{p1,p2})*VB(fEpsG, v2?, p1?)
   *VB(fEpsGStar, v3?, p2?!{p1,p2})*VB(fEpsG, v4?, p2?)
 =  (-d_(v1,v2) - 2*DS(p1-p2,0,1)*(p1(v1)*p2(v2) + p1(v2)*p2(v1)))
   *(-d_(v3,v4) - 2*DS(p1-p2,0,1)*(p1(v3)*p2(v4) + p1(v4)*p2(v3)));
id  VB(fEpsAStar, v1?, p1?{p1,p2})*VB(fEpsA, v2?, p1?)
   *VB(fEpsAStar, v3?, p2?!{p1,p2})*VB(fEpsA, v4?, p2?)
 =  (-d_(v1,v2) - 2*DS(p1-p2,0,1)*(p1(v1)*p2(v2) + p1(v2)*p2(v1)))
   *(-d_(v3,v4) - 2*DS(p1-p2,0,1)*(p1(v3)*p2(v4) + p1(v4)*p2(v3)));

******

b VB;
.sort
Keep Brackets;

* Contract remaining gluon (only after pairs of gluons are contracted)
id VB(fEpsGStar, v1?, p1?)*VB(fEpsG, v2?, p1?) = -d_(v1,v2);
id VB(fEpsAStar, v1?, p1?)*VB(fEpsA, v2?, p1?) = -d_(v1,v2);

******
.sort

#if (`leftNLoops' == 1) && (`dimScheme' == HVIR)
*    b d_,p1,p2,p3,p4,p5,p6,G,e_;
    .sort
    dimension 4;
    Keep Brackets;
    #do i=1,1
*        if ( match(d_(v1?vhvext,v2?$temp)*d_(v2?,v3?)) );
*            sum $temp;
*            redefine i "0";
        if ( match(p1?(v1?vhvext$temp)) );
            sum $temp;
            redefine i "0";
        endif;
*        b d_,p1,p2,p3,p4,p5,p6,G,e_;
        .sort
        Keep Brackets;
    #enddo
    .sort
    dimension d;

#endif * `leftNLoops' == 1 && `dimScheme' == HVIR

b G, G5;
.sort
F g1,g2;
Keep Brackets;

* Gamma5 substitution, up to 3 fermion lines (change hard-coding later)

#do i = 0,`SPINLINES'-1
    repeat id G5(`i')*G(`i') = G5(`i');
    repeat id G(`i', v1?)*G(`i') = G(`i', v1);
    repeat id G5(`i')*G(`i', v1?) = -G(`i',v1)*G5(`i');
    repeat id G5(`i')*G5(`i') = G(`i');

#if (`dimScheme' == HVIR)
    #do j=1,1
        repeat id G(`i',v1?,?a,v2?!{v1?},v1?,?b) = 2*d_(v1,v2)*G(`i',v1,?a,?b)-G(`i',v1,?a,v1,v2,?b);
        repeat id G(`i',?a,v1?,v1?,?b) = d_(v1,v1)*G(`i',?a,?b);
        repeat id G(`i',?a)*G(`i',?b)*G5(`i') = G(`i',?a,?b)*G5(`i');
        id once G(`i',?a,v1?,?b,v1?,?c) = G(`i',?a)*G(`i',v1,?b,v1)*G(`i',?c);
        if ( match(G(`i',v1?,?a,v1?)) ) redefine j "0";
        .sort
    #enddo
    id G(`i',?a)*G5(`i') = distrib_(-1,4,g1,g2,?a);
    id g1(?args) = e_(?args);
    id g2(?args) = g_(`i',?args);
* Any gamma5 left should be alone
    id G5(`i') = 0;
*    repeat id G(`i',v1?,v2?,?args) = G(`i',v1)*G(`i',v2,?args);
#endif * `dimScheme' == HVIR
#enddo

id G(v1?,v2?) = g_(v1,v2);
id G(v1?) = gi_(v1);
id G5(v1?) = g5_(v1);
*id G5(v1?) = 1/24*g_(v1,vdhv0,vdhv1,vdhv2,vdhv3)*e_(vdhv0,vdhv1,vdhv2,vdhv3);
.sort

#if `dimScheme' == HV
* Note on numerical traces: g5 is not yet implemented, and the traces have
* not been implemented in the C code. Do not use for now.
    #if `NUMERICALTRACES' == 1
* Chisholm all but first trace, so remaining traces will have no connected
* indices.
**        #do i = 1, `SPINLINES'-1
**            chisholm `i';
**        #enddo
        b g_,g5_;
        .sort
        cf trace;
        Keep Brackets;
* Apply some simple identities to eliminate repeated indices, similarly
* to how FORM traces internally.
        #do i = 0, `SPINLINES'-1
**            repeat;
**                repeat;
**                    repeat id g_(`i',v2?)*g_(`i',v2?) = 4*gi_(`i');
**                    repeat id g_(`i',p1?)*g_(`i',p1?) = p1.p1*gi_(`i');
**                    id g_(`i',v2?)*g_(`i',v3?)*g_(`i',v4?)*g_(`i',v2?) = 4*gi_(`i')*d_(v3,v4);
**                    id g_(`i',v2?)*g_(`i',v3?)*g_(`i',v4?)*g_(`i',v5?)*g_(`i',v2?) = -2*g_(`i',v5)*g_(`i',v4)*g_(`i',v3);
**                endrepeat;
*** Build a chain with whatever is left, and apply more general identities
**                id once ifnomatch->11 g_(`i',v2?dummyindices_$temp) = G($temp);
**                    repeat id G(?args)*g_(`i',v2?) = G(?args,v2);
**                    repeat id g_(`i',v2?)*G(?args) = G(v2,?args);
**                    if ( multipleof(2) != count(G,1) ) id G(?args) = 0;
**                    repeat id G(v1?,?args,$temp,?args2,$temp,?args3) = g_(`i',v1)*G(?args,$temp,?args2,$temp,?args3);
**                    repeat id G(?args,$temp,?args2,$temp,?args3,v1?) = G(?args,$temp,?args2,$temp,?args3)*g_(`i',v1);
*** While we have the whole chain, check first/last pairs
**                    repeat;
**                        id G(v1?,v2?,?args,v1?) = 4*G(v2,?args);
**                        id G(p1?,v2?,?args,p1?) = p1.p1*G(v2,?args);
**                    endrepeat;
**                    chainout G;
**                    if ( multipleof(2) == count(G,1) );
**                        chainin G;
**                        id G($temp,v1?,?args,v2?,$temp) = 2*G(v2,?args)+2*G(reverse_(?args),v2);
**                    else;
**                        chainin G;
**                        id G($temp,?args,$temp) = -2*G(reverse_(?args));
**                    endif;
**                    chainout G;
**                    id G(v1?) = g_(`i',v1);
**                label 11;
**            endrepeat;
            b g_,g5_;
            .sort
            Keep Brackets;
            mul trace();
            repeat id trace(?args)*g_(`i',v1?) = trace(?args, v1);
            repeat;
                id trace(p1?,v1?,?args,p1?) = p1.p1*trace(v1,?args);
                id trace(?args,v1?,p1?,p1?,?args2) = p1.p1*trace(?args,v1,?args2);
                id trace(?args,p1?,p1?,v1?,?args2) = p1.p1*trace(?args,v1,?args2);
            endrepeat;
            id trace(v1?) = 0;
            id trace()*gi_(`i') = 4;
            id trace() = 1;
            id trace(v1?,v2?) = d_(v1,v2);
        #enddo
* Temporary: for speed testing, delete gammas.
*        id g_(v1?,v2?) = 1;
    #else
        #do i = 0, `SPINLINES'-1
            trace4, `i';
        #enddo
        .sort
    #endif * `NUMERICALTRACES' == 1
    .sort

#else
    #call Traces(3)

*** Why would we need to sum an auto-sum index??
*    #do i = 0,9
*        sum vdhv`i';
*    #enddo
    .sort

    dimension 4;
    #do i = 0,19
        sum vhv`i'l,vhv`i'r;
        sum v`i'el;
    #enddo
    .sort

#endif * `dimScheme' == HV

* Final epsilon tensor manipulation

b e_;
.sort

Keep Brackets;
#message contracting eps tensor
contract, 0;

b e_;
*print;
.sort 
Keep Brackets;

* This must be replaced with the option for the full expression of an
* epsilon tensor (this has been done in the FourVector class in tred).
#if `disableEpsTensor' == 1
    id e_(p1?,p2?,p3?,p4?) = 0;
#endif

b DS;
.sort
Keep Brackets;

repeat id once DS(p1?,m1?,n1?)^n2? = DS(p1,m1,n1*n2);
repeat id once DS(p1?,m1?,n1?)*DS(p1?,m1?,n2?) = DS(p1,m1,n1+n2);

#call eliminateMomentum
#include processSpecific.inc

* New, more general DS elimination
b DS, p1,p2,p3,p4,p5,p6;
.sort
Keep Brackets;

splitarg DS;
id DS(-p1?,p2?,m1?,n1?) = DS(p1,-p2,m1,n1);
repeat id once DS(?args,n1?)^n2? = DS(?args,n1*n2);
repeat id once DS(?args,n1?)*DS(?args,n2?) = DS(?args,n1+n2);

repeat;
*    id once DS(p1?, -p2?, m1?, n1?)*p1?.p2? = -1/2*(DS(p1,-p2,m1,n1-1)-m1^2*DS(p1,-p2,m1,n1));
*    id once DS(p1?, p2?, m1?, n1?)*p1?.p2? = 1/2*(DS(p1,p2,m1,n1-1) +m1^2*DS(p1,p2,m1,n1));
    id once DS(p1?, -p2?, 0, n1?)*p1?.p2? = -1/2*(DS(p1,-p2,0,n1-1)-(mass(p1)^2+mass(p2)^2)*DS(p1,-p2,0,n1));
    id once DS(p1?, p2?, 0, n1?)*p1?.p2? = 1/2*(DS(p1,p2,0,n1-1)-(mass(p1)^2+mass(p2)^2)*DS(p1,p2,0,n1));
    id DS(?args, 0) = 1;
endrepeat;

#include processSpecific.inc
* I think this repeated module was needed in an old version where a dot
* product could be reintroduced after processSpecific. Test if needed and
* eliminate if not during cleanup.

b DS, p1,p2,p3,p4,p5,p6;
.sort
Keep Brackets;

**repeat;
***    id once DS(p1?, -p2?, m1?, n1?)*p1?.p2? = -1/2*(DS(p1,-p2,m1,n1-1)-m1^2*DS(p1,-p2,m1,n1));
***    id once DS(p1?, p2?, m1?, n1?)*p1?.p2? = 1/2*(DS(p1,p2,m1,n1-1) +m1^2*DS(p1,p2,m1,n1));
**    id once DS(p1?, -p2?, 0, n1?)*p1?.p2? = -1/2*(DS(p1,-p2,0,n1-1)-(mass(p1)^2+mass(p2)^2)*DS(p1,-p2,0,n1));
**    id once DS(p1?, p2?, 0, n1?)*p1?.p2? = 1/2*(DS(p1,p2,0,n1-1)-(mass(p1)^2+mass(p2)^2)*DS(p1,p2,0,n1));
**    id DS(?args, 0) = 1;
**endrepeat;

repeat id DS(p1?,p2?,?args) = DS(p1+p2,?args);

.sort

#endprocedure
