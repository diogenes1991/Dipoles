#-
off statistics;
#include declarations.h
#include config.inc
#include `rightDiagDefFile' # global

#do i = `first', `last'
    load `savDir'/simpdiagr`i'.sav;
#enddo

l sumLO = 
#do i = `first', `last'
    + simpdiagr`i'
#enddo
;

* What is this doing here?? How can there be a Lorentz contraction loop at LO?? (Seth)
b d;
.sort
id d = 4 - 2*ep;

b colr, i_, `bracketconstants';
.sort
Keep Brackets;
#write <`idDir'/sumLO.out> "sumLO = %E;", sumLO
.end
