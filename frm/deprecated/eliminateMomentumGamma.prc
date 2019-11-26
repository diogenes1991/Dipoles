#procedure eliminateMomentumGamma
    b VB;
    .sort
    Keep Brackets;
    id VB(n1?, p2?, p1?, ?args) = VB(n1, acc(p2), acc(p1), ?args);
    id VB(n1?, v1?, p1?, ?args) = VB(n1, acc(v1), acc(p1), ?args);

    b p1,p2,p3,p4,p5;
    .sort
    id `ELIMMOM';
    argument;
        id `ELIMMOM';
    endargument;
    b VB;
    .sort
    Keep Brackets;

    id VB(n1?, acc(p2?), acc(p1?), ?args) = VB(n1, p2, p1, ?args);
    id VB(n1?, acc(v1?), acc(p1?), ?args) = VB(n1, v1, p1, ?args);
    b G;
    .sort
    Keep brackets;
    splitarg G;
    repeat id G(n1?, p1?, p2?, ?args) = G(n1,p1) + G(n1,p2,?args);
    id G(n1?,-p1?) = -G(n1,p1);
    .sort
#endprocedure
