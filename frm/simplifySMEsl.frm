#-
off statistics;
#include declarations.h
#include config.inc
#include `leftDiagDefFile' # global

#include `idDir'/SMEsl.expr
b SCC;
.sort
Keep Brackets;
id SCC(?args) = SC(?args);

b SC;
.sort
Keep Brackets;
repeat id once SC(n1?,n2?,?args) = n2*SC(n1,?args);
id SC(n1?) = 1;
b QQ;
.sort
Keep Brackets;
id QQ(?args, l) = QQ(?args);
id QQ(?args, r) = QQ(?args);
.sort
cf tag, indtag;

b G;
.sort
Keep Brackets;
*#call simplifyFiveGammas
*#call simplifyFiveGammas
*#call simplifyFiveGammas
*#call simplifySMEadvanced

#call eliminateMomentum
#call diracAlgebra
#call applyDiracEquation
#call diracAlgebra
#call simplifyBosons

#call collectSC;

mul dummySME;
b SCC, VB, dummySME;
.sort

id dummySME*SCC(?args) = SCC(?args);
id dummySME*VB(?args) = VB(?args);
.sort

b SCC, VB, e_, dummySME;
print;
.end
