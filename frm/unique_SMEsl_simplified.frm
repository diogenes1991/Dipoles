#-
off statistics;
#include declarations.h
#include config.inc
cf SimpSMEl;

#include `idDir'/unique_SMEsl_simplified.expr

b SCC, VB, e_, dummySME;
.sort
cf f;
Keep Brackets;
collect f;

id f(n?) = 1;
.sort

b SCC, VB, e_, dummySME;
print;
.end

