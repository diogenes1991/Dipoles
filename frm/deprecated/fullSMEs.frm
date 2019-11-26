#-
off statistics;
#include declarations.h
#include config.inc
table, sparse, LOpart(1);
table, sparse, SME(2);
cf SMEl, SMEr, SimpSMEl, SimpSMEr;

#include `idDir'/LOparts.id
#include `idDir'/SMEProducts.id

#include `idDir'/fullSMEDefinitions.inc

b SimpSMEl, SimpSMEr;
.sort
Keep Brackets;
id SimpSMEl(n1?)*SimpSMEr(n2?) = SME(n1, n2);

b DS;
.sort
Keep Brackets;
id DS(p1?,m1?,n1?) = DS(p1,m1,1)^n1;
normalize, (0), DS;
*splitarg DS;
*b DS, p1, p2, p3, p4, p5, p6;
*.sort
*Keep Brackets;
*
*repeat;
*    id DS(p1?, p2?, m?, n?)*p1?.p2? = 1/2*DS(p1,p2,m,n-1) - 1/2*(p1.p1 + p2.p2 + m^2)*DS(p1,p2,m,n);
*    id DS(p1?, ?args, m?, 0) = 1;
*endrepeat;
*
*repeat id DS(p1?,p2?,?args,m?, n?) = DS(p1 + p2, ?args,m, n);
b DS;
.sort
Keep Brackets;
repeat id DS(p1?, m?,n1?)*DS(p1?,m?, n2?) = DS(p1,m,n1+n2);

#include processSpecific.inc;


b p1,p2,p3,p4,p5,p6;
.sort
format mathematica;
Keep Brackets;

b p1,p2,p3,p4,p5,p6;
.sort
Keep Brackets;

id p1?.p2? = SS(p1,p2);
*#do i = 1,6
*    #do j = `i',6
*        #do n = 2,5
*            id (p`i'.p`j')^`n' = p`i'p`j'pow`n';
*        #enddo
*        id p`i'.p`j' = p`i'p`j';
*    #enddo
*#enddo

b d;
.sort
Keep Brackets;
id d = 4 - 2*ep;

b ep;
.sort
Keep Brackets;
#do i = 0,{`nSMEl' -1} 
    #do j = 0, {`nLOparts'-1}
        l [fullSME({`i'*`nLOparts' + `j'},0)] = [fullSME({`i'*`nLOparts' + `j'})][1];
        l [fullSME({`i'*`nLOparts' + `j'},1)] = [fullSME({`i'*`nLOparts' + `j'})][ep];
        l [fullSME({`i'*`nLOparts' + `j'},2)] = [fullSME({`i'*`nLOparts' + `j'})][ep^2];
    #enddo
#enddo
.sort
#do i = 0, {`nSMEl' -1 }
    #do j = 0, {`nLOparts' -1}
        hide [fullSME({`i'*`nLOparts' + `j'})];
    #enddo
#enddo

* Change by Seth: factor out masses for more compact SMEs.
b DS, ep, i_, `bracketconstants', `bracketmasses';
print;
.end
