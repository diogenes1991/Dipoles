#-
#:MaxTermSize 100K
off statistics;
#include declarations.h
#include config.inc

#define diagfile1 "`leftDiagDefFile'"
#define diagfile2 "`rightDiagDefFile'"
#define suffix1 "l"
#define suffix2 "r"
.global

#do i = 1,2
    #include `diagfile`i'' # global
    #if `subsetNLODiagrams' == 1 && `i' == 1
        #define first "`subsetfirst'"
        #define last "`subsetlast'"
    #endif
    cf acc;
    cf dtag;
    .global

    #do j = `first', `last'
        load sav/diag`suffix`i''`j'.sav;
    #enddo

* This certainly explains why this script uses so much memory. Apparently,
* the method by which identical structures between diagrams are identified
* is to sum all of the diagrams. This works, but perhaps it would be better
* to let Python identify unique SMEs and color structures. Also, "dtag" is
* never used and is dropped in every case by DropCoefficient.
    #message summing up
    g expr`i' = 
    #do j = `first', `last'
        + dtag(`j')*diag`suffix`i''`j'
    #enddo
    ;
    .store;

    #message color
    l dd = expr`i';
    b color, delta, cOlOne, delta3;
    .sort
    collect acc,acc,40;
    DropCoefficient;
    id acc(?args) = 1;
    DropCoefficient;
    b color;
    .sort
    Keep Brackets;
    id color(n1?) = n1;
*    b cOlT, cOlf, cc, delta, cOlOne, delta3;
** CR 20170330: Also bracket cOlTr, for new diagram-level color treatment.
    b cOlT, cOlf, cc, delta, cOlOne, delta3, cOlTr;
** CR
    .sort
    Keep Brackets;
    #write<`idDir'/color`suffix`i''.out> "color = %E;", dd;
    .store

    #message spinor chains

    l dd = expr`i';
    b SCC, VB, e_;
    .sort
    collect acc,acc,40;
    DropCoefficient;
    id acc(?args) = 1;
    DropCoefficient;

    b SCC, VB, e_;
    .sort
    Keep Brackets;
    #write<`idDir'/SMEs`suffix`i''.out> "smes = %E;", dd;

    b SCC;
    .sort
    collect acc,acc,40;
    DropCoefficient;
    id acc(?args) = 1;
    DropCoefficient;

    b SCC;
    .sort
    Keep Brackets;
    #write<`idDir'/SCCs`suffix`i''.out> "sccs = %E;", dd;
    #close<`idDir'/SCCs`suffix`i''.out>
    .store

#enddo
.end
