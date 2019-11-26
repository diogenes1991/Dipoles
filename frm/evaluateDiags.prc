#procedure evaluateDiags(FILE, MARKER, NUM, NUMPIECES, CONJ)

dimension d;
#include `FILE' # global
.global
#if `subsetNLODiagrams' == 1 && `LOOPS' == 1
    #define first "`subsetfirst'"
    #define last "`subsetlast'"
#endif
* Make sure we don't redo diagrams
#if `NUM' >= `NUMPIECES'
    #redefine last "0"
#endif



#do i = `first'+`NUM', `last', `NUMPIECES'
    #include `FILE' # d`i'
    .store

    g diag`MARKER'`i' = d`i';
    .sort:test message;
    #message d`i'

    #call treatDiagram(`CONJ')
    #call collectSC
    b VB, SCC;
    .sort:First collectSC;
    Keep Brackets;

*   * Strip color
    id VB(n1?, v1?, p1?, cOli1?) = VB(n1,v1,p1);
    repeat id SCC(?args1, QQ(?args2, p1?, cOli1?, n2?), ?args4) = 
              SCC(?args1, QQ(?args2, p1, n2), ?args4);

    #include `kinematicsFile'
*    #call colorMatrices
** CR 20170330: Alternatively there is another routine for handling color on the
** diagram level.  May be switched through the input-file switch useOldDiagColor
** (true/false) corresponding to oldDiagColor (1/0).
    #if `oldDiagColor' == 1
    #call colorMatrices
    #elseif `oldDiagColor' == 0
    #call NLOXSimplifyColorMatrices
    #endif
** CR

* It is now safe to absorb cOlOne

    id cOlOne^n? = cOlOne;
    id cOlOne*cOlT(cOli1?,cOli2?,cOli3?) = cOlT(cOli1,cOli2,cOli3);
    id cOlOne*cOlf(cOli1?,cOli2?,cOli3?) = cOlf(cOli1,cOli2,cOli3);
    id cOlOne*delta(cOli1?,cOli2?) = delta(cOli1,cOli2);

    .sort:Color matrices;

** Conjugate diagram couplings/masses if in complex mass scheme
*    #if `CONJ' == 1 && `complexMass' == 1
*        id m1?partmasses?conjmasses = m1;
*        id n1?ccouplings?ccouplingsconj = n1;
    #if `CONJ' == 1
        repeat id m1?partmasses?conjmasses = m1;
        repeat id n1?ccouplings?ccouplingsconj = n1;
        repeat id (m1?partmasses?conjmasses)^(-1) = (m1)^(-1);
        repeat id (n1?ccouplings?ccouplingsconj)^(-1) = (n1)^(-1);
    #endif

*   * This #if block only operates on closed fermion loops, and the #if statement
*   * should be modified to reflect this. 
    #if `LOOPS' == 1
*      * CR, 20170211: Comment out the following code block. For the case of
*      * maximally one fermion loop, being the case at one-loop, this is too
*      * complicated, and it also breaks processes with no external spin lines
*      * but fermion loops. Instead, we can simply identify each "loose" G and
*      * G5 with a spin-line index of 0 then to be traced over. Each G and G5
*      * on an external spin line is caught in an SCC and won't be affected by
*      * this until released later on by which time the trace over the spin-line
*      * index 0 has been already taken.
*        b SCC;
*        .sort:Bracket SCC;
*        Keep Brackets;
*        id SCC(n1?, ?args) = acc(n1)*SCC(n1,?args);
*
*        b G, G5, acc;
*        .sort;Bracket Gs;
*        Keep Brackets;
*
*       * Changed to just one pass (can only have one fermion loop at NLO)
*        #do i=0,1
*            id acc(n1?)*G(n2?!{n1?},v1?) = acc(n2,`i')*G(n2, v1);
*            id acc(n1?)*G(n2?!{n1?}) = acc(n2,`i')*G(n2);
*            id acc(n1?)*G5(n2?!{n1?}) = acc(n2,`i')*G5(n2);
*        #enddo
*        id acc(n1?) = 1;
*        repeat id acc(n1?, n3?)*acc(n1?,n3?) = acc(n1, n3);
*        b acc,G,G5 ;
*        .sort:Test message;
*        Keep Brackets;
*        repeat id acc(n1?,n2?)*G(n1?,v1?) = acc(n1,n2)*G(n2,v1);
*        repeat id acc(n1?,n2?)*G(n1?) = acc(n1,n2)*G(n2);
*        repeat id acc(n1?,n2?)*G5(n1?) = acc(n1,n2)*G5(n2);
*        id acc(n1?,n2?) = 1;

*       * CR, 20170211: As stated above, the following will replace above's
*       * code block. With these two lines we still get the correct numbers
*       * for uubar_ddbar_QCD and at the same time gg_gg with fermion loops
*       * compiles.

        id G(n1?,?args) = G(0,?args);
        id G5(n1?) = G5(0);
*        #call diracAlgebra

        b G, G5;
        .sort:Dirac algebra;
        F g1,g2,test;
        Keep Brackets;

*       * #do loop unnecessary; only i = 0 possible
        #do i = 0,1
            repeat id G5(`i')*G(`i') = G5(`i');
            repeat id G(`i', v1?)*G(`i') = G(`i', v1);
            repeat id G5(`i')*G(`i', v1?) = -G(`i',v1)*G5(`i');
            repeat id G5(`i')*G5(`i') = G(`i');
	    #if `dimScheme' == HV || `dimScheme' == HVIR
                #do j=1,1
                    repeat id G(`i',v1?,?a,v2?!{v1?},v1?,?b) = 
                              2*d_(v1,v2)*G(`i',v1,?a,?b)-G(`i',v1,?a,v1,v2,?b);
                    repeat id G(`i',?a,v1?,v1?,?b) = d_(v1,v1)*G(`i',?a,?b);
                    repeat id G(`i',?a)*G(`i',?b)*G5(`i') = G(`i',?a,?b)*G5(`i');
                    id once G(`i',?a,v1?,?b,v1?,?c) =
                            G(`i',?a)*G(`i',v1,?b,v1)*G(`i',?c);
                    if ( match(G(`i',v1?,?a,v1?)) ) redefine j "0";
                    .sort:G5 loop;
                #enddo

                id G(`i',?a)*G5(`i') = distrib_(-1,4,g1,g2,?a);
                id g1(?args) = e_(?args);
                id g2(?args) = g_(`i',?args);

*               * Any gamma5 left should be alone
                id G5(`i') = 0;
*               * repeat id G(`i',v1?,v2?,?args) = G(`i',v1)*G(`i',v2,?args);
	    #endif
        #enddo

        id G(v1?,v2?) = g_(v1,v2);
        id G(v1?) = gi_(v1);
        id G5(v1?) = g5_(v1);

        b g_, e_;
        .sort:Gamma algebra;

        #if `dimScheme' == HV || `dimScheme' == HVIR
*           * changed by Seth to d-dim trace
            #call Traces(1);

        #endif * `dimScheme' == HV

    #endif * `LOOPS' == 1

*   * Unpack the SCC objects
    b SCC;
    .sort:Trace, SCC brackets;
    Keep Brackets;
    id SCC(?args) = SC(?args);

    b SC;
    .sort:SC bracket;
    on properorder;
    Keep Brackets;
    repeat id once SC(n1?,n2?,?args) = n2*SC(n1,?args);
    id SC(n1?) = 1;
    b QQ;
    .sort:SC simplification;
    Keep Brackets;
    id QQ(?args, l) = QQ(?args);
    id QQ(?args, r) = QQ(?args);
    .sort:QQ l/r removal;

    b e_;
    .sort:Epsilon handling;
    cf tag, indtag;
    Keep Brackets;
    contract, 0;

    b G;
    .sort:Epsilon contraction;
    b VB;
    .sort:VB bracket;
    Keep Brackets;

    id VB(n1?, p1?, p2?, ?args) = VB(n1, acc(p1), p2, ?args);
    .sort:VB momentum preservation;

    b VB;
    .sort:VB index relabeling;
    Keep Brackets;
    id VB(n1?, acc(p1?), p2?, ?args) = VB(n1, p1, p2, ?args);

    b G;
    .sort:G bracket;
    Keep Brackets;

    #call eliminateMomentum
    #call diracAlgebra
    #call simplifyBosons
    #call applyDiracEquation
    #call diracAlgebra
    #call simplifyBosons

    #if `LOOPS' == 1 && (`dimScheme' == HVIR)

*       * Relabel 4d indices
        b d4d;
        .sort
        Keep Brackets;
        repeat;
            id d4d(p1?,v1?) = d_(p1,v1);
** The commented out statements in this block should be redundant now that
** d4d is declared symmetric. These changes require further testing.
*            id d4d(v1?,p1?) = d_(p1,v1);
            id d4d(v1?,v2?)*d4d(v1?,v3?) = d4d(v2,v3);
*            id d4d(v1?,v2?)*d4d(v2?,v3?) = d4d(v1,v3);
*            id d4d(v1?,v2?)*d4d(v3?,v1?) = d4d(v2,v3);
*            id d4d(v1?,v2?)*d4d(v3?,v2?) = d4d(v1,v3);
        endrepeat;
        .sort
        Keep Brackets;
        #do i=0,2
            .sort
            Keep Brackets;
            if ( match(d4d(v1?,v2?!ve$temp)) );
                multiply replace_($temp,v`i'el);
                id d4d(v1?,v`i'el) = d_(v1,v`i'el);
            endif;
        #enddo
    #endif

    id G(n1?)*G(n1?,v1?) = G(n1,v1);
    id G(n1?)*G5(n1?) = G5(n1);
    id G(n1?,v1?)*G(n1?) = G(n1,v1);
    id G5(n1?)*G(n1?) = G5(n1);
*   *  #call simplifyBosons
    id QQ(?args, l) = QQ(?args);
    id QQ(?args, r) = QQ(?args);

    #call eliminateMomentum
    #call simplifyBosons
    #include `kinematicsFile';
    .sort:Bosons and kinematics;

*    ab cc, cOlT, cOlf, delta, cOlOne, delta3;
** CR 20170330: Also bracket cOlTr, for new diagram-level color treatment.
    ab cc, cOlT, cOlf, delta, cOlOne, delta3, cOlTr;
** CR
    .sort:color antibracket;
    collect color;
    splitarg color;
    repeat id color(n1?, n2?, ?args) = color(n1) + color(n2, ?args);
    b color;
    .sort:color bracket;
    Keep Brackets;
    factarg, (-1), color;
    id color(?args,n2?, n1?) = n1*color(?args, n2);

*   * sDS is a denominator created by treatDiagrams.prc. From what I can
*   * tell, it stores any denominator which doesn't contain a loop-
*   * momentum. Couldn't the replacement of sDS to DS be done there?
    b sDS;
    .sort:sDS bracket;
    Keep Brackets;

    id sDS(?args) = DS(?args);
    b DS;
    .sort:DS bracket;
    Keep Brackets;
*   Note: here we conjugate all denominator masses. We only conjugate
*   numerator masses in the complex mass scheme; check if this is correct.
    #if `CONJ'
        id DS(p1?, m1?partmasses?conjmasses, n1?) = DS(p1, m1, n1);
    #endif
    repeat id DS(p1?, m1?, n1?)*DS(p1?,m1?,n2?) = DS(p1,m1,n1+n2);

    #call collectSC
    b SCC;
    .sort:SCC bracket;
    cf tags;
    Keep Brackets;

*   * Use smallest momentum index carrying fermion to sort chains 
    id SCC(n1?, QQ(?args3,p1?,?args), ?args2) =
       SCC(n1, QQ(?args3,p1,?args), ?args2)*tags(p1);
    id tags(p1?)*tags(p2?) = tags(p1, p2);
    repeat id tags(p1?, ?args1)*tags(p2?, ?args2) = tags(p1, ?args1,p2,?args2);
    symmetrize tags;
    b SCC, tags;
    .sort:SCC tags;
    Keep Brackets;

* Changed by Seth: add third Dirac chain handling
*    id SCC(n1?, QQ(?args3, p1?,?args), ?args2)*tags(p1?, ?args3) =
*        SCC(n0`MARKER', QQ(?args3,p1,?args), ?args2)*tags(p1, ?args3);
*    id SCC(n1?, QQ(?args3, p2?,?args), ?args2)*tags(p1?, p2?, ?args3) =
*        SCC(n1`MARKER', QQ(?args3,p2,?args), ?args2)*tags(p1, p2, ?args3);
*    id SCC(n1?, QQ(?args3, p3?,?args), ?args2)*tags(p1?, p2?, p3?) = 
*        SCC(n2`MARKER', QQ(?args3,p3,?args), ?args2);
** CR 20170216: The above six lines break processes with more than 2 spin lines.
** Using ?args3 twice on the left side in the first two id statements above will
** actually cause them to not be triggered. This is a bug: in each one of those,
** the first ?args3 and the second ?args3 are supposed to be different arguments
** and we need to use different argument list wildcards. The purpose of the con-
** struct is as follows: if by any chance, in the same set of diagrams some dia-
** grams would be assigned a different ordering in their spin-line assigment wrt
** the particles involved, then this should re-assign those values consistently.
** NOTE though, that the whole thing only works for up to 3 external spin lines.
    id SCC(n1?, QQ(?args3, p1?,?args), ?args2)*tags(p1?, ?args4) =
        SCC(n0`MARKER', QQ(?args3,p1,?args), ?args2)*tags(p1, ?args4);
    id SCC(n1?, QQ(?args3, p2?,?args), ?args2)*tags(p1?, p2?, ?args4) =
        SCC(n1`MARKER', QQ(?args3,p2,?args), ?args2)*tags(p1, p2, ?args4);
    id SCC(n1?, QQ(?args3, p3?,?args), ?args2)*tags(p1?, p2?, p3?) = 
        SCC(n2`MARKER', QQ(?args3,p3,?args), ?args2);

*    id SCC(n1?, QQ(n2?, p1?,?args), ?args2)*tags(p1?, ?args3) =
*        SCC(n0`MARKER', QQ(n2,p1,?args), ?args2)*tags(p1, ?args3);
*    id SCC(n1?, QQ(n2?, p2?,?args), ?args2)*tags(p1?, p2?) =
*        SCC(n1`MARKER', QQ(n2,p2,?args), ?args2);

*   if there is only one dirac chain, remove tags()
    id tags(?args) = 1;
*    id tags(p1?) = 1;

    repeat id SCC(n1?, ?args1, G(n2?!{,n1?}, v1?), ?args2) =
        SCC(n1, ?args1, G(n1,v1,0), ?args2);
    repeat id SCC(n1?, ?args1, G5(n2?!{,n1?}), ?args2) =
        SCC(n1, ?args1, G5(n1,0), ?args2);
    repeat id SCC(n1?, ?args1, G(n2?!{,n1?}), ?args2) =
        SCC(n1, ?args1, G(n1,0), ?args2);
    repeat id SCC(n1?, ?args1, PL(n2?!{,n1?}), ?args2) =
        SCC(n1, ?args1, PL(n1,0), ?args2);
    repeat id SCC(n1?, ?args1, PR(n2?!{,n1?}), ?args2) =
        SCC(n1, ?args1, PR(n1,0), ?args2);

    argument;
        id G(n1?,v1?,0) = G(n1,v1);
        id G5(n1?,0) = G5(n1);
        id G(n1?,0) = G(n1);
        id PL(n1?,0) = PL(n1);
        id PR(n1?,0) = PR(n1);
    endargument;

    b SCC, tags;
    .sort:Tagging done;

*    b delta, cOlOne, delta3, i_, VB, SCC, DS, color, TI, p1, p2, p3, p4, p5, `bracketconstants';
    b delta, cOlOne, delta3, i_, VB, SCC, DS, color, TI, p1, p2, p3, p4, p5, p6, `bracketconstants';
    .sort:Final bracket;
    Keep Brackets;
    #write<math/diag`MARKER'`i'.m> "(%E)" diag`MARKER'`i';
    #close<math/diag`MARKER'`i'.m>
*    b delta, cOlOne, delta3, e, i_, g, VB, SCC, DS, color, TI, p1, p2, p3, p4, p5;
    b delta, cOlOne, delta3, e, i_, g, VB, SCC, DS, color, TI, p1, p2, p3, p4, p5, p6;
    .store

    save sav/diag`MARKER'`i'.sav diag`MARKER'`i';
    delete storage; * cleanup before store file gets too big to fit in memory
#enddo * Main do loop

#endprocedure
