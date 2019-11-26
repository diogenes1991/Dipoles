#procedure collectSC

* assign spinor chain to gamma5
b G, G5;
.sort
Keep Brackets;
id G5 = G5(-1);

id G5(n1?)*G(n2?, v1?) = G5(n2)*G(n2,v1);
id G(n2?,v1?)*G5(n1?) = G(n2,v1)*G5(n2);

.sort
s l, r;
s args1, args2, args3, args4;

* In below, all quarks set requirements have been removed in QQ

b QQ;
.sort
Keep Brackets;
#do i=1,10
  id once QQ(?args, p1?, cOli1?) = QQ(?args, p1, cOli1, (-1)^(`i'),0);
  id once QQ(?args, p1?) = QQ(?args, p1, (-1)^(`i'),0);
#enddo

id QQ(?args,p1?, cOli1?,-1,0) = QQ(?args,p1, cOli1, l);
id QQ(?args,p1?, cOli1?,1,0) = QQ(?args,p1, cOli1, r);
id QQ(?args,p1?, -1,0) = QQ(?args,p1, l);
id QQ(?args,p1?, 1,0) = QQ(?args,p1, r);

ab QQ;
.sort:QQ labeling;
#if `VERSION_' >= 4
on oldfactarg;
#endif
collect SC;
splitarg SC;
repeat id once SC(n1?, ?args1, n2?) = SC(n1) + SC(?args1, n2);
factarg,  SC;
id SC(?args, n1?!{,-1}, n2?) = SC(?args, n1, 1, n2);
id SC(?args, n1?, n2?) = n1*n2*SC(?args);
b SC;
.sort:SC isolation;
Keep Brackets;

repeat id SC(QQ(?args1), ?args2, QQ(?args3,r), QQ(?args4,l), ?args5) =
   SC(QQ(?args1), ?args2, QQ(?args3,r))*SC(QQ(?args4,l), ?args5);
    
repeat;
    id SC(?args, G(?args2)) = SC(?args)*G(?args2);
    id SC(?args, G5(?args2)) = SC(?args)*G5(?args2);
endrepeat;
* Change to eliminate constant SC (fermion-less diagrams)
id SC(n1?) = n1;
repeat id once SC(?args1, G(n1?,?args), ?args2) = SCC(n1, ?args1, G(n1,?args), ?args2);
repeat id once SC(?args1, G5(n1?), ?args2) = SCC(n1, ?args1, G5(n1), ?args2);
repeat id once SC(?args1, PL(n1?), ?args2) = SCC(n1, ?args1, PL(n1), ?args2);
repeat id once SC(?args1, PR(n1?), ?args2) = SCC(n1, ?args1, PR(n1), ?args2);
b SCC;
.sort:SC labeling;
Keep Brackets;

id SCC(n1?, QQ(?args1),QQ(?args2)) = SCC(n1, QQ(?args1), G(n1), QQ(?args2));
.sort:SC unit gamma insert;
b SCC;
.sort:SCC bracket;


#endprocedure
