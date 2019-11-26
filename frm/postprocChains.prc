#procedure postprocChains

* here we assume that all spinor chain indices of the right hand side
* have been eliminated!!!

b G, G5;
.sort
Keep Brackets;
#message post processing chains

#do i=0,10
    id G(n`i'l,v1?) = Gc(`i', v1);
    id G(n`i'l) = Gc(`i');
    id G5(n`i'l) = G5c(`i');
#enddo

id Gc(n?, -k?) = -Gc(n, k);

b Gc, G5c;
.sort
Keep Brackets;

id G5c(n1?) = G5(n1);
id Gc(n1?) = G(n1);
id Gc(n1?, v1?) = G(n1, v1);
#endprocedure
