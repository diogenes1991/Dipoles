#procedure preprocChains

b SCC;
.sort
cf dd;
cf SCCd;
Keep Brackets;
#message Preproc - stage 1

repeat id once SCC(n1?, ?args, QQ(n3?, p1?, r))*SCC(n2?, QQ(Bar(n3?), p1?,l), ?args4) =
  SCCd(dd(n1,n2), ?args, QQ(n3, p1,r), QQ(Bar(n3), p1,l), ?args4);

b SCCd;
.sort
Keep Brackets;

#message Preproc - stage 1
repeat id once SCCd(n1?, ?args, QQ(n3?, p1?,?args2, r))*SCCd(n2?, QQ(Bar(n3?), p1?,?args3,l), ?args4) =
  SCCd(n1,n2, ?args, QQ(n3, p1,?args2,r), QQ(Bar(n3), p1,?args3,l), ?args4);

id SCCd(?args) = SCC(?args);

b SCC;
.sort
Keep Brackets;
#message rotating
id SCC(n3?, QQ(n1?, p1?, ?args1), ?args2, QQ(n2?, p1?, ?args3)) = SCC(n3, ?args2, QQ(n2, p1, ?args3), QQ(n1, p1, ?args1));
id SCC(n3?,n4?, QQ(n1?, p1?, ?args1), ?args2, QQ(n2?, p1?, ?args3)) = SCC(n3, n4, ?args2, QQ(n2, p1, ?args3), QQ(n1, p1, ?args1));
b SCC;
.sort

#endprocedure
