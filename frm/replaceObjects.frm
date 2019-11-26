#-
off statistics;
#include declarations.h
#include config.inc

#include `leftDiagDefFile' # global
cf acc;
cf SMEl, SMEr, SimpSMEl, SimpSMEr;
cf colr, coll;
.global

#if `subsetNLODiagrams' == 1
    #define first "`subsetfirst'"
    #define last "`subsetlast'"
#endif

#ifndef `diagNum'
#define diagNum "0"
#endif

#define first "`first'+`diagNum'"

#do i = `first', `last', `numPieces'
    load `savDir'/diag`leftMarker'`i'.sav;
#enddo
.sort

#do i = `first', `last', `numPieces'
  #message diag`i'
  g simpdiag`leftMarker'`i' = + diag`leftMarker'`i';

  b SCC, VB, e_;
  .sort
  Keep Brackets;
  #include `idDir'/SMEs`leftMarker'.id
  #include `idDir'/SMEs`leftMarker'ToSimp.id
  .sort

  b SCC, VB, SimpSMEl;
  .sort
  b color;
  .sort
  Keep Brackets;
  id color(n1?) = n1;

*  b cOlT, cOlf, delta, cOlOne, cc, delta3;
** CR 20170330: Also bracket cOlTr, for new diagram-level color treatment.
  b cOlT, cOlf, delta, cOlOne, cc, delta3, cOlTr;
** CR
  .sort
  Keep Brackets;
  #include `idDir'/color`leftMarker'.id

  b SMEl, SMEr, SimpSMEr, SimpSMEl, e_, coll, colr, DS, TI, delta3, delta, cOlOne, `bracketconstants';
  .sort
  Keep Brackets;

  #write<`mathDir'/simpdiag`leftMarker'`i'.m> "(%E)" simpdiag`leftMarker'`i';
  #close<`mathDir'/simpdiag`leftMarker'`i'.m>
  .store

  save `savDir'/simpdiag`leftMarker'`i'.sav simpdiag`leftMarker'`i';
#enddo
.clear

#-
off statistics;
#include declarations.h
#include config.inc
#include `rightDiagDefFile' # global
cf acc;
cf SMEl, SMEr, SimpSMEl, SimpSMEr;
cf colr, coll;
.global

#do i = `first', `last'
    load `savDir'/diag`rightMarker'`i'.sav;
#enddo
.sort

#do i = `first', `last'
  #message diag`rightMarker'`i'
  g simpdiag`rightMarker'`i' = + diag`rightMarker'`i';

  b SCC, VB;
  .sort
  Keep Brackets;
  #include `idDir'/SMEs`rightMarker'.id
  #include `idDir'/SMEs`rightMarker'ToSimp.id

  b color;
  .sort
  Keep Brackets;
  id color(n1?) = n1;
*  b cOlT, cOlf, delta, cOlOne, cc, delta3;
** CR 20170330: Also bracket cOlTr, for new diagram-level color treatment.
  b cOlT, cOlf, delta, cOlOne, cc, delta3, cOlTr;
** CR
  .sort
  Keep Brackets;
  #include `idDir'/color`rightMarker'.id

  b SMEl, SimpSMEr, SimpSMEl, coll, colr, DS, TI, delta, cOlOne, delta3, `bracketconstants';
  .sort
  Keep Brackets;

  #write<`mathDir'/simpdiag`rightMarker'`i'.m> "(%E)" simpdiag`rightMarker'`i';
  #close<`mathDir'/simpdiag`rightMarker'`i'.m>
  .store

  save `savDir'/simpdiag`rightMarker'`i'.sav simpdiag`rightMarker'`i';
#enddo
.clear
.end
