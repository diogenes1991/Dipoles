#procedure simplifyFiveGammas
b G;
.sort
Keep Brackets;
    id G(n1?, v1?)*G(n1?,v2?)*G(n1?,v3?)*G(n1?,v4?)*G(n1?,v5?) =
       +d_(v1,v2)*G(n1,v3)*G(n1,v4)*G(n1,v5)
       -d_(v1,v3)*G(n1,v2)*G(n1,v4)*G(n1,v5)
       +d_(v1,v4)*G(n1,v2)*G(n1,v3)*G(n1,v5)
       -d_(v1,v5)*G(n1,v2)*G(n1,v3)*G(n1,v4)
       +d_(v2,v3)*G(n1,v1)*G(n1,v4)*G(n1,v5)
       -d_(v2,v4)*G(n1,v1)*G(n1,v3)*G(n1,v5)
       +d_(v2,v5)*G(n1,v1)*G(n1,v3)*G(n1,v4)
       +d_(v3,v4)*G(n1,v1)*G(n1,v2)*G(n1,v5)
       -d_(v3,v5)*G(n1,v1)*G(n1,v2)*G(n1,v4)
       +d_(v4,v5)*G(n1,v1)*G(n1,v2)*G(n1,v3)
       -(d_(v1,v2)*d_(v3,v4) - d_(v1,v3)*d_(v2,v4) + d_(v1,v4)*d_(v2,v3))*G(n1, v5)
       +(d_(v1,v2)*d_(v3,v5) - d_(v1,v3)*d_(v2,v5) + d_(v1,v5)*d_(v2,v3))*G(n1, v4)
       -(d_(v1,v2)*d_(v4,v5) - d_(v1,v4)*d_(v2,v5) + d_(v1,v5)*d_(v2,v4))*G(n1, v3)
       +(d_(v1,v3)*d_(v4,v5) - d_(v1,v4)*d_(v3,v5) + d_(v1,v5)*d_(v3,v4))*G(n1, v2)
       -(d_(v2,v3)*d_(v4,v5) - d_(v2,v4)*d_(v3,v5) + d_(v2,v5)*d_(v3,v4))*G(n1, v1);
#endprocedure
