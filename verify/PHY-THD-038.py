# O20163-1 Otto cycle vs cycle with full expansion to initial pressure (monatomic)
import sympy as sp, numpy as np
from scipy.integrate import quad
r,al,p0,V0=sp.symbols('r alpha p0 V0',positive=True)
g=sp.Rational(5,3)
VA=V0; pA=p0; VB=V0/r; pB=p0*r**g; pC=al*pB; VC=VB
pD=pC*(VC/VA)**g
VF=VC*(pC/pA)**(1/g)
print("pD/pA=",sp.simplify(pD/pA),"  VF/VA=",sp.simplify(VF/VA))
Qin=sp.Rational(3,2)*VB*(pC-pB)
QoutOtto=sp.Rational(3,2)*VA*(pD-pA)
QoutC=sp.Rational(5,2)*pA*(VF-VA)
ea=sp.simplify(1-QoutOtto/Qin); ec=sp.simplify(1-QoutC/Qin)
print("ea=",ea); print("ec=",ec)
print("ec target check:",sp.simplify(ec-(1-sp.Rational(5,3)*r**sp.Rational(-2,3)*(al**sp.Rational(3,5)-1)/(al-1))))
# numeric with integration of p dV for cycle c
def work_c(rv,av):
    p0n=1;V0n=1;VBn=V0n/rv;pBn=p0n*rv**(5/3);pCn=av*pBn;VFn=VBn*(pCn/p0n)**(3/5)
    Wcomp=quad(lambda V:p0n*(V0n/V)**(5/3),V0n,VBn)[0]
    Wexp=quad(lambda V:pCn*(VBn/V)**(5/3),VBn,VFn)[0]
    Wiso=p0n*(V0n-VFn)
    Q=1.5*VBn*(pCn-pBn)
    return (Wcomp+Wexp+Wiso)/Q
for rv,av in [(4,2),(8,3),(2,1.5)]:
    print(rv,av,work_c(rv,av), float(ec.subs({r:rv,al:av})), float(ea.subs({r:rv})))
# subcycle A->D->F->A work positive: W = Q_AD - Q_FA
W_sub=sp.simplify(QoutOtto-QoutC)
print("Wsub=",W_sub)
x=sp.symbols('x',positive=True)
f=x-1-sp.Rational(5,3)*(x**sp.Rational(3,5)-1)
print("f(1)=",f.subs(x,1)," f'(1)=",sp.diff(f,x).subs(x,1)," f''=",sp.simplify(sp.diff(f,x,2)))
print("limit alpha->1 of ratio:",sp.limit((al**sp.Rational(3,5)-1)/(al-1),al,1))
# numbers r=8, alpha=32  -> alpha^(3/5)=8
print("r=8,alpha=32: ea=",ea.subs(r,8)," ec=",sp.nsimplify(ec.subs({r:8,al:32})), " VF/VA=",sp.simplify((VF/VA).subs(al,32)))
