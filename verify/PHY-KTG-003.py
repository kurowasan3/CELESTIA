# P30-1 平均自由行程と熱伝導（連続領域とクヌーセン領域）
import sympy as sp
k,T,p,d,vb,L,dT=sp.symbols('k T p d vbar L DeltaT',positive=True)
nden=p/(k*T)
lam=sp.solve(sp.Eq(nden*sp.pi*d**2*sp.Symbol('lam'),1),sp.Symbol('lam'))[0]; print('lambda=',lam)
f=sp.Rational(1,4)*nden*vb
# model A: molecules come from distance lam on each side, carry 3/2 k T
grad=dT/L
qA=sp.simplify(f*sp.Rational(3,2)*k*(2*lam*grad)); print('qA=',qA,' p-dependence:',sp.diff(qA,p))
qB=sp.simplify(f*sp.Rational(3,2)*k*dT); print('qB=',qB)
ratio=sp.simplify(qB/qA); print('qB/qA=',sp.simplify(ratio), ' check L/(2lam):',sp.simplify(ratio-L/(2*lam)))
pc=sp.solve(sp.Eq(qA,qB),p)[0]; print('crossover p=',pc)
# numeric: lambda0=6.8e-8 at 1.0e5 Pa; L=5.0e-3; want q = qA/100 in B regime -> L/(2 lam)=1/100
lam0,p0,Lv=6.8e-8,1.0e5,5.0e-3
lam_needed=50*Lv; pneed=p0*lam0/lam_needed; print('p needed=',pneed)
# lambda at crossover numeric
print('p_c numeric =',p0*lam0/(Lv/2))
# check 1/100 by direct formulas with Ar-like numbers
import math
kB=1.38e-23;Tv=300; dv=math.sqrt(kB*Tv/(math.pi*p0*lam0)); vbv=400
qAv=float(qA.subs({k:kB,T:Tv,d:dv,vb:vbv,L:Lv,dT:1}))
qBv=float(qB.subs({k:kB,T:Tv,p:pneed,vb:vbv,dT:1}))
print('d=',dv,'qB/qA=',qBv/qAv)
# flux relative to 1 atm at 1.0e2 Pa and 1.0 Pa (actual flux = min(qA,qB))
for pv in [1.0e2,1.0]:
    lamv=lam0*p0/pv
    rel=min(1.0, Lv/(2*lamv))
    print('p=',pv,'lambda=',lamv,'q/q_atm=',rel)
