#procedure diracAlgebraIndexed(ind)

b G, GG;
.sort
dimension d;
Keep Brackets;
id G(`ind')*G(`ind',v2?)  = G(`ind',v2);
id G(`ind',v2?)*G(`ind')  = G(`ind',v2);
repeat;
    id disorder G(`ind'?,v1?)*G(`ind'?,v2?) = 2*d_(v1,v2) - G(`ind',v2)*G(`ind',v1);
    id disorder G(`ind'?,p1?)*G(`ind'?,p2?) = 2*p1.p2 - G(`ind',p2)*G(`ind',p1);
*    id GG(5)*G(`ind'?,v1?) = -G(`ind',v1)*GG(5);
*    id GG(5)*GG(5) = 1;
    id G(`ind'?,v2?)*G(`ind'?,v2?) = 4;
    id G(`ind'?,p1?)*G(`ind'?,p1?) = p1.p1;
endrepeat;
.sort

#endprocedure
