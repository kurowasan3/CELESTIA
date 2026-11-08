# S0206-1 verification: two-cylinder regenerative engine, rectangle cycle in (V_H, V_C)
import sympy as sp
import numpy as np
from scipy.integrate import quad
R,T,V0=sp.symbols('R T V0',positive=True)
VH,VC=sp.symbols('V_H V_C',positive=True)
TH,TC=3*T,T
# 1 mol total: p*VH/(R TH)+p*VC/(R TC)=1
p=sp.solve(sp.Eq(sp.Symbol('p')*VH/(R*TH)+sp.Symbol('p')*VC/(R*TC),1),sp.Symbol('p'))[0]
print('p =',sp.simplify(p))
# steps: 1(V0,V0)->2(4V0,V0)->3(4V0,2V0)->4(V0,2V0)->1
Y12=sp.integrate(p.subs(VC,V0),(VH,V0,4*V0)); Y34=sp.integrate(p.subs(VC,2*V0),(VH,4*V0,V0))
Z23=sp.integrate(p.subs(VH,4*V0),(VC,V0,2*V0)); Z41=sp.integrate(p.subs(VH,V0),(VC,2*V0,V0))
Y=sp.simplify(sp.logcombine(Y12+Y34,force=True)); Z=sp.simplify(sp.logcombine(Z23+Z41,force=True))
print('Y =',Y,' expect',3*R*T*sp.log(sp.Rational(49,40)), sp.simplify(Y-3*R*T*sp.log(sp.Rational(49,40))))
print('Z =',Z,' Z/Y =',sp.simplify(Z/Y))
# heat from hot body via dQ = p dVH - R TH dnH
nH=sp.simplify(p*VH/(R*TH))
states=[(V0,V0),(4*V0,V0),(4*V0,2*V0),(V0,2*V0),(V0,V0)]
nHs=[sp.simplify(nH.subs({VH:a,VC:b})) for a,b in states]
print('nH at states',nHs)
Qs=[Y12-R*TH*(nHs[1]-nHs[0]), -R*TH*(nHs[2]-nHs[1]), Y34-R*TH*(nHs[3]-nHs[2]), -R*TH*(nHs[4]-nHs[3])]
print('step Q_H:',[sp.nsimplify(sp.simplify(q/(R*T))) for q in Qs])
QH=sp.simplify(sum(Qs)); print('QH-Y =',sp.simplify(QH-Y))
eta=sp.simplify((Y+Z)/QH); print('eta =',eta)
# entropy check: Q_C (released) / TC == Q_H / TH
QC=-Z; print('QH/TH - QC/TC =',sp.simplify(QH/TH-QC/TC))
# general rectangle a->a', b->b'
a,a2,b,b2=sp.symbols('a a2 b b2',positive=True)
Yg=sp.integrate(p.subs(VC,b),(VH,a,a2))+sp.integrate(p.subs(VC,b2),(VH,a2,a))
Zg=sp.integrate(p.subs(VH,a2),(VC,b,b2))+sp.integrate(p.subs(VH,a),(VC,b2,b))
print('general Z/Y =',sp.simplify(sp.expand_log(Zg,force=True)/sp.expand_log(Yg,force=True)))
# numeric check of the path integral of heat (fine steps) with numbers
Rn,Tn,V0n=8.31,300.0,1e-3
pf=lambda h,c: 3*Rn*Tn/(h+3*c)
def seg(fh,fc,N=200000):
    s=np.linspace(0,1,N+1); h=fh(s); c=fc(s); pm=pf((h[1:]+h[:-1])/2,(c[1:]+c[:-1])/2)
    nh=pf(h,c)*h/(Rn*3*Tn)
    return np.sum(pm*np.diff(h))-Rn*3*Tn*(nh[-1]-nh[0]), np.sum(pm*np.diff(h)), np.sum(pm*np.diff(c))
segs=[(lambda s:V0n*(1+3*s),lambda s:V0n+0*s),(lambda s:4*V0n+0*s,lambda s:V0n*(1+s)),(lambda s:V0n*(4-3*s),lambda s:2*V0n+0*s),(lambda s:V0n+0*s,lambda s:V0n*(2-s))]
res=np.array([seg(*f) for f in segs]); QHn,Yn,Zn=res.sum(axis=0)
print('numeric QH,Y,Z',QHn,Yn,Zn,' eta',(Yn+Zn)/QHn,' formula Y',3*Rn*Tn*np.log(49/40))
