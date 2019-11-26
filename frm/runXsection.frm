#-
#:SmallSize 50M
#:LargeSize 200M
#:ScratchSize 500M
#:MaxTermSize 100K
off statistics;
#include declarations.h
#include config.inc

s Pi;

#include `leftDiagDefFile' # global
#ifndef `diagNum'
#define diagNum "0"
#endif
#define firstl "`first'+`diagNum'"
#define lastl  "`last'"
#if `subsetNLODiagrams' == 1
    #define firstl "`subsetfirst'+`diagNum'"
    #define lastl "`subsetlast'"
#endif

cf SMEl, SMEr, SimpSMEl, SimpSMEr;
cf coll, colr;
*cf SMESum;
table,sparse, SME(2);
table,sparse, fullSME(1);
table,sparse, col(2);
table,sparse, LOpart(2);
#include `idDir'/fullSMEsNull.id
#include `idDir'/ColorProducts.id
.global

#include `idDir'/LOdef.id
.store
dimension d;

#do i = `firstl', `lastl', `numPieces'
    #message d`i'
    load sav/simpdiagl`i'.sav;
    .store
    #message mult2
    g dd`i' =  simpdiagl`i'*LOdef*`xsectNorm';

* Initial-state color average
    mul cOlNR^(-`nInQuarks');
    mul cOlNA^(-`nInGluons');

* Initial-state spin average
    mul 2^(-`nInFermions');
    mul 2^(-`nInMasslessVBs');
    mul 3^(-`nInMassiveVBs');

* Final-state symmetry factor
    mul 1 / `symmetryFactor';

* Reduce powers of sqrt2
    id sqrt2^2 = 2;
    id sqrt2^-2 = 1/2;

*   added by Seth; may need to change for d dimensions
*(1+`nInBosons'*ep+`nInBosons'*(`nInBosons'+1)/2*ep^2);

*   added by Seth; multiply by renormalization factor if on
    #if `renorm' == 1
        mul `renormFactor';
    #endif

    b colr, coll;
    .sort
    Keep Brackets;
    #message Colors
    id coll(n1?)*colr(n2?) = col(n1,n2);

    b cOlNR, cOlNA, cOlcA, cOlcR, cOlOne;
    .sort
    Keep Brackets;
    id cOlcR^n? = CF^n;
    id cOlcA^n? = CA^n;
    id cOlNR^n? = NF^n;
    id cOlNA^n? = NA^n;
    id cOlOne = 1;

    b CF, CA, NF, NA, color;
    .sort 
    Keep Brackets;

    argument;
* Temporary to fix missing color factor (Steve: wtf is this?)
*        id cOld44(cOlpA1?cOlpAs,cOlpA2?cOlpAs) = 1/48*NF^2*(-72 + 70*NF^2 + 2*NF^4);
        id CF^n1? = (4/3)^n1;
        id CA^n1? = 3^n1;
        id NA^n1? = 8^n1;
        id NF^n1? = 3^n1;
    endargument;
    id CF^n1? = (4/3)^n1;
    id CA^n1? = 3^n1;
    id NA^n1? = 8^n1;
    id NF^n1? = 3^n1;

    id color(n1?)*CF^n2? = color(n1*CF^n2);
    id color(n1?)*CA^n2? = color(n1*CA^n2);
    id color(n1?)*NA^n2? = color(n1*NA^n2);
    id color(n1?)*NF^n2? = color(n1*NF^n2);
    id color(n1?) = n1;

    b SMEl, SimpSMEl, LOpart;
    .sort 
    Keep Brackets;
    #message SMEs
    id SimpSMEl(n1?)*LOpart(n2?) = fullSME(n1*`nLOparts' + n2);
    id LOpart(n2?) = fullSME(n2);
*    id SMEl(n1?)*LOpart(n2?) = sum_(n,0,2,fullSME(n1,n2,n)*ep^n);

    b TI;
    .sort
    Keep Brackets;

*    id TI(?args) = sum_(n,-2,0,TI(?args,n)*ep^n);

    b DS, p1, p2, p3, p4, p5, p6;
    .sort
    Keep Brackets;

* Changed by Seth: probable incorrect sign on m term below
    repeat;
        id DS(p1?, p2?, m?, n?)*p1?.p2? = 1/2*DS(p1,p2,m,n-1) - 1/2*(p1.p1 + p2.p2 - m^2)*DS(p1,p2,m,n);
        id DS(p1?, ?args, m?, 0) = 1;
    endrepeat;

    repeat id DS(p1?,p2?,?args,m?, n?) = DS(p1 + p2, ?args,m, n);
    repeat id DS(p1?, m?,n1?)*DS(p1?,m?, n2?) = DS(p1,m,n1+n2);

    #include processSpecific.inc;

    b p1,p2,p3,p4,p5,p6;
    .sort
    Keep Brackets;
    id p1?.p2? = SS(p1,p2);
*
*    #do i = 1,6
*        #do j = `i',6
*            #do n = 2,5
*                id (p`i'.p`j')^`n' = p`i'p`j'pow`n';
*            #enddo
*            id p`i'.p`j' = p`i'p`j';
*        #enddo
*    #enddo

*   * Why would there still be sDS objects at this stage? They
*   * should have been eliminated in evaluateDiags.prc.
    b sDS;
    .sort
    Keep Brackets;
    id sDS(?args) = DS(?args);
    b DS;
    .sort 
    Keep Brackets;
    id DS(p1?,m1?,n1?) = DS(p1,m1,1)^n1;
    b DS;
    .sort
    Keep Brackets;
    normalize, (0), DS;
    repeat id DS(p1?, m?, n1?)*DS(p1?,m?, n2?) = DS(p1,m,n1+n2);
    #include processSpecific.inc

    b d;
    .sort
    format mathematica;
    Keep Brackets;
    id d = 4-2*ep;

    b ep;
* remove print to keep log file smaller
*    print;
    .sort

*   Collect all ep dependence
    ab ep;
    .sort
    collect acc;
    normalize acc;
    b acc;
* remove print to keep log file smaller
*    print;
    .sort
    Keep Brackets;
    splitarg acc;
    repeat id acc(n1?, n2?, ?args) = acc(n1)*acc(n2, ?args);
    mul acc(1)*acc(ep)*acc(ep^2);
    factarg,(-1),acc;
    b acc;
    .sort
    Keep Brackets;

    id acc(ep,n1?)*acc(ep,n2?) = acc(ep,n1*n2)^2;
*   Added by Seth, guess to keep ep^2 terms?
*    id acc(ep^2,n1?)*acc(ep^2,n2?) = acc(ep^2,n1*n2)^2;
    mul epc(0,0,0);
    b acc, epc;
    .sort
    Keep Brackets;
    id acc(1)^2*epc(0,0,0) = epc(1,0,0);
    id acc(ep,n?)^2*epc(n1?,n2?,n3?) = epc(n1,n,n3);
    id acc(ep^2,n?)^2*epc(n1?,n2?,n3?) = epc(n1,n2,n);
    id acc(?args) = 1;

*   Get rid of constants expressed as Poly3
    id epc(1,0,0) = 1;

*   Collect SMEs to keep expressions smaller (by Seth)
    ab fullSME;
    .sort

    #message SME sum extraction
    collect ISMESum,ISMESum,40;

    MakeInteger ISMESum;

    b ISMESum;
    .sort

*   minimize necessary ep-dependent multiplications in final numerical
*   code by pulling ep dependent coefficient out of the whole
*   expression. (Temp disabled by Seth)
*    b DS, TI, color, `bracketconstants', `bracketcolor', epc;
*    .sort
*    collect acc;
*    normalize acc;

    Keep Brackets;

    #$sumNum = -1;
    $sumNum = $sumNum+1;
    id ISMESum(?args) = ISMESum($sumNum,?args);
    #message ismesum replace
    .sort

    L stripDD = dd`i';
    b ISMESum;
    .sort

* shorten expression a bit so next collect succeeds
*    skip dd`i';
*    dropcoefficient;
*    id epc(?args) = 1;
*    b SMESum;
*    .sort

    skip dd`i';
    collect acc;
    id acc(?args) = 1;
    .sort

    skip dd`i';
    DropCoefficient;
    .sort

    #write <`idDir'/ISMEsums`i'.out> "ismesums = %E;", stripDD;
    #close <`idDir'/ISMEsums`i'.out>
    .sort

    b ISMESum;
    .sort

    id ISMESum(n?, ?args) = ISMESum(n);

* Group dot product expressions
* Note: in principle it would probably be cleaner to put the mass logs
* in with the prefactors (constants), but this would require rewriting
* the pickling of the mass log variables in genCCode and testing.
* In any case they only appear in counterterms, which are not the
* large/slow part.
*    ab SS, `bracketmasses', `bracketmasslogs';
    ab SS, `bracketmasses';
    .sort

    #message Grouping Dot Products
    collect SSSum,SSSum,40;

    MakeInteger SSSum;

    b SSSum;
    .sort

* Remove this next line if we want to handle this in genCCode
    id SSSum(1) = 1;

    ab `bracketconstants';
    #message constcollect
    .sort
#if `VERSION_' >= 4
    on OldFactArg;
#endif
*    cf const;
    collect const;
    normalize const;
    splitarg const;
    repeat id const(n1?,n2?,?args) = const(n1) + const(n2,?args);
    normalize const;
    factarg const;
    repeat id const(n1?,n2?,?args) = const(n1)*const(n2, ?args);
    id const(1) = 1;
    repeat;
        id const(n1?)*const(n1?) = const(n1^2);
        id const(n1?)*const(n1?^n3?) = const(n1^(1+n3));
        id const(n1?^n2?)*const(n1?^n3?) = const(n1^(n2+n3));
    endrepeat;
    repeat;
        id const(n1?)*const(n1?) = const(n1^2);
        id const(n1?)*const(n1?^n3?) = const(n1^(1+n3));
        id const(n1?^n2?)*const(n1?^n3?) = const(n1^(n2+n3));
    endrepeat;

* Collect prefactors with common expressions into one
    ab DS, color, `bracketconstants', `bracketcolor', const;
    .sort
    collect dum_;
    normalize dum_;

    b dum_, DS, color, `bracketconstants', `bracketcolor', const, TI, pole;
    print +s;
    .sort
    Keep Brackets;
    #write <`mathDir'/xsectiond`i'.m> "[xsectiond`i'] = %E;\n\n", dd`i';
    #close <`mathDir'/xsectiond`i'.m>
    .store
    save `savDir'/xsectiond`i'.sav dd`i';
#enddo

.end
