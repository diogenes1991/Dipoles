#procedure simplifyColor
.sort
#message >> simplify color
dimension cOlNA;
b delta, cOlOne, delta3, cOlT, cOlf;
.sort
Keep Brackets;
id delta(cOli1?, cOli2?) = d_(cOli1, cOli2);
id cOlOne = 1;
repeat;
    id delta3(cOli1?, cOli2?)*delta3(cOli2?, cOli3?) = delta3(cOli1, cOli3);
* temp by Seth, seem to be missing a case
*    id delta3(cOli1?, cOli2?)*delta3(cOli3?, cOli1?) = delta3(cOli2, cOli3);
    id delta3(cOli1?, cOli1?) = 3;
    id delta3(cOli2?, cOli1?)*cOlT(cOli1?, cOli3?, cOli4?) = cOlT(cOli2, cOli3, cOli4);
    id delta3(cOli2?, cOli1?)*cOlT(cOli3?, cOli1?, cOli4?) = cOlT(cOli3, cOli2, cOli4);
endrepeat;
* Temp comment by Seth to see if there are delta3s left at this point
*id delta3(cOli1?, cOli2?) = d_(cOli1, cOli2);

b cOlT, cOlf;
.sort
Keep Brackets;

#call docolor
.sort
dimension d;

id cOlI2R^n? = (cOlNR/cOlNA*CF)^n;

id cOlcR^n? = CF^n;
id cOlcA^n? = CA^n;
id cOlNR^n? = NF^n;
id cOlNA^n? = NA^n;

*id CF*NF/NA = 1/2;

id TF/CF = NF/NA;

id cOld33(cOlpR1?cOlpRs,cOlpR2?cOlpRs) = CF/8*(NF^2 - 4);
id cOld44(cOlpR1?cOlpRs,cOlpR2?cOlpRs) = CF/48*NF^-1*(18 - 6*NF^2 + NF^4);
id cOld44(cOlpR1?cOlpRs,cOlpA1?cOlpAs) = 1/48*NF*(-6 + 5*NF^2 + NF^4);
id cOld44(cOlpA1?cOlpAs,cOlpA2?cOlpAs) = 1/48*NF^2*(-72 + 70*NF^2 + 2*NF^4);
id cOld55(cOlpR1?cOlpRs,cOlpR2?cOlpRs) = CF/384*NF^-2*(NF^2 - 2)*(NF^4 + 24);
id cOld433(cOlpR1?cOlpRs,cOlpR2?cOlpRs,cOlpR3?cOlpRs) = CF/96*NF^-1*(NF^2 - 4)*(NF^2 - 6);

ab nf, CF, CA, TF, NF, NA,
   cOlpA1, ..., cOlpA`RANK',
   cOlpR1, ..., cOlpR`RANK',
   cOld33,cOld44,cOld55,cOld433,cOld66,cOld633,cOld543,cOld444,cOld3333,
   cOld77,cOld743,cOld653,cOld644,cOld554,cOld5333,cOld4433a,cOld4433b,cOld4433c,
   cOld88,cOld853,cOld844,cOld763,cOld754,cOld7333,cOld664,cOld655,cOld6433a,cOld6433b,cOld6433c,
   cOld5533a,cOld5533b,cOld5533c,cOld5533d,cOld5443a,cOld5443b,cOld5443c,cOld5443d,
   cOld4444a,cOld4444b,cOld43333a,cOld43333b,
   cOlf, cOlT;
.sort

collect color;
normalize color;
.sort

#endprocedure
