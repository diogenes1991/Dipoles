**************************************************
*                                                *
* declaration of common symbols                  *
*                                                *
**************************************************

* general
**********

s d, ep;
dimension d;
.sort
cf acc;

auto i v;
auto i vhv=0;

#do i=0,19
  i v`i'el=0;
  i v`i'er=0;
#enddo

auto s p1p,p2p,p3p,p4p,p5p,p6p;

set ve:<v0el>,...,<v19el>,<v0er>,...,<v19er>;
* set vhv:vhv0,...,vhv19;
set vhvext:<vhv0l>,...,<vhv19l>,<vhv0r>,...,<vhv19r>,<v0el>,...,<v19el>;

* vdhv are indices connecting two different spin lines in the same left diagram,
* to be preserved as d-dimensional
set vdhv:vdhv0,...,vdhv19;

* Kronecker for color
cf delta(Symmetric);
cf delta3(Symmetric);
s cOlOne; * Used for colorless fields

* Temporary d_ tensor for 4d (by Seth)
ct d4d(symmetric);

* momenta (up to 6-point functions)
* CR -- not sure whether this gets in conflict with auto s p
* CR -- why is no autodeclare used here
v P, p1, p2, p3, p4, p5, p6;
v k;
* v K,K1, K2, K3, K4, K5, K6;
v k1, k2, k3, k4, k5, k6;

set extmom:p1,p2,p3,p4,p5,p6;

* masses
s M, m, m1, m2, m3, m4, m5, m6;
* TODO: For the particle-specific masses, like many of the declarations in
* this file, we really ought to build these automatically from the model.
s mu, md, mb, mt, mW, mZ, mH;
s must, mdst, mbst, mtst, mWst, mZst, mHst;
set partmasses:mu,md,mb,mt,mW,mZ,mH;
set conjmasses:must,mdst,mbst,mtst,mWst,mZst,mHst;

s args;

s pole;


* propagators
f SF;
cf DS, DG, DV;
cf SS(Symmetric);
cf SSinv(Symmetric);

* vertices
*
* VVV(k1,a1,k2,a2,k3,a3) = (k1(a3) - k2(a3))*d_(a1, a2) + ...
*
* S(k) = k.k
*
* Hold for loop momenta
* KK(k(mu)) = k(mu);
*
* T(k, a1, a2) = k.k*d_(a1,a2) - k(a1)*k(a2)
*
* G(i,a) = g_(i,a)

cf VVV, VVO;
cf S, KK;
f G, G5, G4;

* gauge parameter (Feynman gauge: xi = 0)
* Possibly something that works? Untested
s xi(:1);

* flavours and SU(3) invariants
s nlf, nhf;
s CF, CA, TF, NF, NA;
cf color;

* color package (Vermaseren) 
#include color.h

.sort
dimension cOlNA;
s cOlNA, cOlNR;
#do i = 0, 30
    i cOli`i'e;
    i cOli`i'l;
    i cOli`i'r;
    i cOlj`i'l;
    i cOlj`i'r;
#enddo
*auto i cOli=cOlNA;
*auto i cOlj=cOlNA;
set cOlie:<cOli0e>,...,<cOli29e>;
.sort
dimension d;

s N;
auto s n;
set nl:<n0l>,...,<n9l>;
set nr:<n0r>,...,<n9r>;

cf sDS, TT;

f SC, SSP;
cf SCC;
cf cc;

** for the color correlators
** try with hard coded symbols first
** also include below in #define bracketconstants and set ccouplings
s ccL0,ccL1,ccL2,ccL3,ccL4,ccL5,ccL6;

f lm, rm;
f func, Gm, ff, ff1, ff2, ff3, ff4;


* external fields

f Ub, UbBar, Vb, VbBar, Uu, UuBar, Vu, VuBar, Ud, UdBar, Vd, VdBar, Ut, VtBar, UqBar, Vq;
f Ul, Vl, VlBar, UlBar;
cf EpsG, EpsGStar, EpsW, EpsWStar, EpsZ, EpsZStar, EpsA, EpsAStar;

s U,V;
s fEpsG, fEpsGStar, fEpsW, fEpsWStar, fEpsZ, fEpsZStar, fEpsA, fEpsAStar;

set gluons: fEpsG, fEpsGStar;
* set bosons: fEpsW, fEpsWStar, fEpsZ, fEpsZStar, fEpsA, fEpsAStar;

set vectorbosons: fEpsG, fEpsGStar, fEpsW, fEpsWStar, fEpsZ, fEpsZStar, fEpsA, fEpsAStar;

f QQ;
f Bar;
cf VB;

* For processes with no external fields/spinor chains
s dummySME;

** CR -- why is no autodeclare used here as well?
v q1,...,q8;

s QWWZ;
s gWH, gZH;
* CR -- not sure whether this gets in conflict with auto s p
cf pol, ds, sme;
cf tt(symmetric);
s Zv, Za;

s Qq, Qt, Qb, Qu, Qd, Ql, e, g;
f Vt, UtBar, VqBar, Uq, EpsAStar;

s l, r;

cf mass;

cf cdelta(Symmetric);

s Qqi, QQi;
s mQ, mq;
* CR -- Already have an auto s n above
*s nf, nif, nff;

*s Qqivv, Qqiav, Qqvv, Qqav, Quvv, Quav, Qdvv, Qdav, Qtvv, Qtav, Qbvv, Qbav, QQiav, QQivv;
*s WudC, WubC, WusC, WtbC;
* new constants for EW; couplings/constants should be auto-added in future
*s g3H,gHchi,gchiHA,gchiHZ,gphiA,gphiZ,gphiHW,gphichiW,gphiA,gphiZ,gphiHW,gphichiW,gphichiW,gWZphi,gWAphi,ychit,ychib,yphivv,yphiav;

* renorm constants
s dZAb, dZmb, dZVb, dZAZ, dZ1Z, dZAb0, dZmb0, dZVb0, dZAZ0, dZ1Z0, dZextZ;
s dZxgsxQ, dZxGxQ, dZxqxQ, dZxtxQ, dZxmtxQ, dZxbxQ, dZxmbxQ;
s dRxGxQ, dRxqxQ, dRxqmxQ, dRxbxQ, dRxtxQ;
s dZxgsxQ0, dZxGxQ0, dZxqxQ0, dZxtxQ0, dZxmtxQ0, dZxbxQ0, dZxmbxQ0;
s dRxGxQ0, dRxqxQ0, dRxqmxQ0, dRxbxQ0, dRxtxQ0;
s dZxgsxQ1, dZxGxQ1, dZxqxQ1, dZxtxQ1, dZxmtxQ1, dZxbxQ1, dZxmbxQ1;
s dRxGxQ1, dRxqxQ1, dRxqmxQ1, dRxbxQ1, dRxtxQ1;

* t'Hooft-Veltman conversion factor
s HV1;

* EW renormalization constants, added by Steve
s dZfLiiE, dZfRiiE, dmfiE;
s dZfLlE, dZfRlE;
s dZfLnE, dZfRnE;
s dZfLuuE, dZfRuuE;
s dZfLddE, dZfRddE;
s dZfLbbE, dZfRbbE, dmbE;
s dZfLttE, dZfRttE, dmtE;
s dZAAE, dZAZE, dZZAE, dZZZE, dZWE, dZHE;
s dmH2E, dmW2E, dmZ2E, dZeE, dswE, dtE;
s dcwE;

s dZfLiiE0, dZfRiiE0, dmfiE0;
s dZfLlE0, dZfRlE0;
s dZfLnE0, dZfRnE0;
s dZfLuuE0, dZfRuuE0;
s dZfLddE0, dZfRddE0;
s dZfLbbE0, dZfRbbE0, dmbE0;
s dZfLttE0, dZfRttE0, dmtE0;
s dZAAE0, dZAZE0, dZZAE0, dZZZE0, dZWE0, dZHE0;
s dmH2E0, dmW2E0, dmZ2E0, dZeE0, dswE0, dtE0;

s dZfLiiE1, dZfRiiE1, dmfiE1;
s dZfLlE1, dZfRlE1;
s dZfLnE1, dZfRnE1;
s dZfLuuE1, dZfRuuE1;
s dZfLddE1, dZfRddE1;
s dZfLbbE1, dZfRbbE1, dmbE1;
s dZfLttE1, dZfRttE1, dmtE1;
s dZAAE1, dZAZE1, dZZAE1, dZZZE1, dZWE1, dZHE1;
s dmH2E1, dmW2E1, dmZ2E1, dZeE1, dswE1, dtE1;

* CR: Need conjugates of fermion dZ's in CM scheme
s dZfLiiEst, dZfRiiEst;
s dZfLlEst, dZfRlEst;
s dZfLnEst, dZfRnEst;
s dZfLuuEst, dZfRuuEst;
s dZfLddEst, dZfRddEst;
s dZfLbbEst, dZfRbbEst;
s dZfLttEst, dZfRttEst;
*
s dZfLiiE0st, dZfRiiE0st;
s dZfLlE0st, dZfRlE0st;
s dZfLnE0st, dZfRnE0st;
s dZfLuuE0st, dZfRuuE0st;
s dZfLddE0st, dZfRddE0st;
s dZfLbbE0st, dZfRbbE0st;
s dZfLttE0st, dZfRttE0st;
*
s dZfLiiE1st, dZfRiiE1st;
s dZfLlE1st, dZfRlE1st;
s dZfLnE1st, dZfRnE1st;
s dZfLuuE1st, dZfRuuE1st;
s dZfLddE1st, dZfRddE1st;
s dZfLbbE1st, dZfRbbE1st;
s dZfLttE1st, dZfRttE1st;

s vu, au, vd, ad, vl, al, vn, an, sw, sw2, cw, cw2, yt, yb, gfp, gfm, sqrt2;
s vust, aust, vdst, adst, vlst, alst, vnst, anst, swst, sw2st, cwst, cw2st, ytst, ybst;
s gup, gum, gdp, gdm, glp, glm, gnm;

set ccouplings: vu,au,vd,ad,vl,al,vn,an,sw,sw2,cw,cw2,yt,yb,ccL0,ccL1,ccL2,ccL3,ccL4,ccL5,ccL6;
set ccouplingsconj: vust,aust,vdst,adst,vlst,alst,vnst,anst,swst,sw2st,cwst,cw2st,ytst,ybst;

s dgupE, dgumE, dgdpE, dgdmE, dglpE, dglmE, dgnmE;
s T3u, T3d, T3l, T3n;

* Things to pull out of overall expressions for simplification.
* In principle we could use sets rather than a preprocessor definition if
* a new enough version of FORM supporting sets in brackets is used.
* Unfortunately to be safe we need to repeat long lists for older versions.
#define bracketconstants "e, g, Qt, Qb, Qd, Qu, Ql, Qq, Qqi, QQi, nlf, nhf, nif, nff, 
dZAb,dZmb,dZVb,dZAZ,dZ1Z,dZAb0,dZmb0,dZVb0,dZAZ0,dZ1Z0,dZextZ,
vu,au,vd,ad,vl,al,vn,an,sw,sw2,cw,cw2,yt,yb,
vust,aust,vdst,adst,vlst,alst,vnst,anst,swst,sw2st,cwst,cw2st,ytst,ybst,
HV1,
dZfLiiE0, dZfRiiE0, dmfiE0,
dZfLlE0, dZfRlE0,
dZfLnE0, dZfRnE0,
dZfLuuE0, dZfRuuE0,
dZfLddE0, dZfRddE0,
dZfLbbE0, dZfRbbE0, dmbE0,
dZfLttE0, dZfRttE0, dmtE0,
dZAAE0, dZAZE0, dZZAE0, dZZZE0, dZWE0, dZHE0,
dmH2E0, dmW2E0, dmZ2E0, dZeE0, dswE0, dtE0,
dZfLiiE1, dZfRiiE1, dmfiE1,
dZfLlE1, dZfRlE1,
dZfLnE1, dZfRnE1,
dZfLuuE1, dZfRuuE1,
dZfLddE1, dZfRddE1,
dZfLbbE1, dZfRbbE1, dmbE1,
dZfLttE1, dZfRttE1, dmtE1,
dZAAE1, dZAZE1, dZZAE1, dZZZE1, dZWE1, dZHE1,
dmH2E1, dmW2E1, dmZ2E1, dZeE1, dswE1, dtE1,
gfp, gfm, sqrt2,
gup, gum, gdp, gdm, glp, glm, gnm,
T3u, T3d, T3l, T3n,
dZfLiiE0st, dZfRiiE0st,
dZfLlE0st, dZfRlE0st,
dZfLnE0st, dZfRnE0st,
dZfLuuE0st, dZfRuuE0st,
dZfLddE0st, dZfRddE0st,
dZfLbbE0st, dZfRbbE0st,
dZfLttE0st, dZfRttE0st,
dZfLiiE1st, dZfRiiE1st,
dZfLlE1st, dZfRlE1st,
dZfLnE1st, dZfRnE1st,
dZfLuuE1st, dZfRuuE1st,
dZfLddE1st, dZfRddE1st,
dZfLbbE1st, dZfRbbE1st,
dZfLttE1st, dZfRttE1st,
dZxgsxQ0, dZxGxQ0, dZxqxQ0, dZxtxQ0, dZxmtxQ0, dZxbxQ0, dZxmbxQ0,
dRxGxQ0, dRxqxQ0, dRxqmxQ0, dRxbxQ0, dRxtxQ0,
dZxgsxQ1, dZxGxQ1, dZxqxQ1, dZxtxQ1, dZxmtxQ1, dZxbxQ1, dZxmbxQ1,
dRxGxQ1, dRxqxQ1, dRxqmxQ1, dRxbxQ1, dRxtxQ1,
ccL0,ccL1,ccL2,ccL3,ccL4,ccL5,ccL6"

#define bracketcolor "NF, NA, CF, CA, cOlI2R, cOlNA"
#define bracketmasses "mu,md,mb,mt,mW,mZ,mH,must,mdst,mbst,mtst,mWst,mZst,mHst"

* Functions for final cross section output
cf const, SSSum, ISMESum, epc;
* Tensor integral coefficient (last in this file so it appears last in output,
* respecting bracket ordering))
cf TI;

f PL, PR;

.global
