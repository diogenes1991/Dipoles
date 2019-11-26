#-
off statistics;
#include declarations.h
#include config.inc

#include `idDir'/unique_SMEsr_simplified.expr

b SCC, VB, dummySME;
.sort
cf f;
Keep Brackets;
collect f;

id f(n?) = 1;
.sort

b SCC, VB, dummySME;
print;
.end
