#procedure eliminateMomentum
    b VB;
    .sort:VB;
    Keep Brackets;
    id VB(n1?, v1?, p1?, ?args) = VB(n1, v1, acc(p1), ?args);

    b p1,p2,p3,p4,p5;
    .sort:momentum bracket;
    id `ELIMMOM';
    argument;
        id `ELIMMOM';
    endargument;
    b VB;
    .sort:VB;
    Keep Brackets;

    id VB(n1?, v1?, acc(p1?), ?args) = VB(n1, v1, p1, ?args);
    b G;
    .sort:G;
    Keep brackets;
    splitarg G;
    repeat id G(n1?, p1?, p2?, ?args) = G(n1,p1) + G(n1,p2,?args);
    id G(n1?,-p1?) = -G(n1,p1);
    b VB;
    .sort:VB;
    Keep Brackets;
    id VB(n1?,v1?, p1?) = VB(n1,v1,acc(p1));
    id VB(n1?,v1?, p1?, cOli1?) = VB(n1,v1,acc(p1), cOli1);
    splitarg VB;

    repeat id once VB(fEpsA,p1?, p2?,?args, acc(p3?)) = VB(fEpsA, p1, acc(p3)) + VB(fEpsA, p2, ?args, acc(p3));
    repeat id once VB(fEpsAStar,p1?, p2?,?args, acc(p3?)) = VB(fEpsAStar, p1, acc(p3)) + VB(fEpsAStar, p2, ?args, acc(p3));

    repeat id once VB(fEpsZ,p1?, p2?,?args, acc(p3?)) = VB(fEpsZ, p1, acc(p3)) + VB(fEpsZ, p2, ?args, acc(p3));
    repeat id once VB(fEpsZStar,p1?, p2?,?args, acc(p3?)) = VB(fEpsZStar, p1, acc(p3)) + VB(fEpsZStar, p2, ?args, acc(p3));

    repeat id once VB(fEpsW,p1?, p2?,?args, acc(p3?)) = VB(fEpsW, p1, acc(p3)) + VB(fEpsW, p2, ?args, acc(p3));
    repeat id once VB(fEpsWStar,p1?, p2?,?args, acc(p3?)) = VB(fEpsWStar, p1, acc(p3)) + VB(fEpsWStar, p2, ?args, acc(p3));

    repeat id once VB(fEpsG,p1?, p2?,?args, acc(p3?)) = VB(fEpsG, p1, acc(p3)) + VB(fEpsG, p2, ?args, acc(p3));
    repeat id once VB(fEpsGStar,p1?, p2?,?args, acc(p3?)) = VB(fEpsGStar, p1, acc(p3)) + VB(fEpsGStar, p2, ?args, acc(p3));

    id VB(n1?,-p1?,acc(p2?)) = -VB(n1,p1,acc(p2));
    id VB(n1?,-p1?,acc(p2?),cOli1?) = -VB(n1,p1,acc(p2),cOli1);
    id VB(n1?,v1?,acc(p2?)) = VB(n1,v1,p2);
    id VB(n1?,v1?,acc(p2?),cOli1?) = VB(n1,v1,p2,cOli1);
#endprocedure
