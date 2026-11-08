# N60-1 muonic atom: Bohr model inside uniformly charged nucleus
import sympy as sp, numpy as np
k,Z,e,R,r,m,hb,n,me,a0,EH,mu=sp.symbols('k Z e R r m hbar n m_e a0 E_H mu',positive=True)
# potential energy inside: integrate force from r to R plus U(R)
F=k*Z*e**2*r/R**3
s=sp.Symbol('s',positive=True)
U_in=-k*Z*e**2/R - sp.integrate((k*Z*e**2*s/R**3),(s,r,R))
print('U_in =',sp.factor(sp.simplify(U_in)))
v=sp.Symbol('v',positive=True)
sol=sp.solve([sp.Eq(m*v**2/r,F), sp.Eq(m*v*r,n*hb)],[v,r],dict=True)[0]
rn=sol[r]; vn=sol[v]; print('r_n =',rn)
En=sp.simplify((m*vn**2/2+U_in).subs(r,rn)); print('E_n =',sp.simplify(En))
w=sp.sqrt(k*Z*e**2/(m*R**3)); print('E_n - (n hbar w - 3kZe^2/2R) =',sp.simplify(En-(n*hb*w-3*k*Z*e**2/(2*R))))
# express with a0, EH: k e^2 = 2 EH a0 ; hbar^2 = me k e^2 a0 ; m = mu me
subs={k:2*EH*a0/e**2}
hbsq=me*2*EH*a0*a0
rn4=sp.simplify((rn**4).subs(hb**2,hbsq).subs(subs).subs(m,mu*me).subs(hb,sp.sqrt(hbsq)))
print('r_n^4 =',sp.simplify(rn4))
hw=sp.simplify((hb*w).subs(subs).subs(m,mu*me).subs(hb,sp.sqrt(hbsq)))
print('hbar w =',hw)
# numbers
num={a0:5.3e-11,EH:13.6,mu:207,Z:82,R:7.1e-15}
r1in=float((rn4.subs(n,1).subs(num))**sp.Rational(1,4)); r2in=float((rn4.subs(n,2).subs(num))**sp.Rational(1,4))
r1pt=5.3e-11/(207*82); r2pt=4*r1pt
print('inside r1,r2 =',r1in,r2in,' point r1,r2=',r1pt,r2pt,' R=',7.1e-15)
hwv=float(hw.subs(num)); E1=hwv-3*13.6*5.3e-11*82/7.1e-15; E2=-207*82**2*13.6/4
E1pt=-207*82**2*13.6
print('hbar w [MeV]',hwv/1e6,' E1in',E1/1e6,' E2out',E2/1e6,' photon',(E2-E1)/1e6,' point photon',(E2-E1pt)/1e6, 'ratio',(E2-E1)/(E2-E1pt))
# independent numeric check in SI
kk=8.99e9; ee=1.602e-19; mme=9.109e-31; hbar=1.0546e-34
mmu=207*mme; Rn=7.1e-15; Zz=82
w_=np.sqrt(kk*Zz*ee**2/(mmu*Rn**3)); print('SI hbar w MeV',hbar*w_/ee/1e6, ' r1',np.sqrt(hbar/(mmu*w_)))
# continuity check of energy at r=R between models? E_out(r)= -kZe^2/(2r); E_in(r)= m w^2 r^2 - 3kZe^2/(2R); at r=R both -kZe^2/(2R)
print(sp.simplify((m*w**2*R**2-3*k*Z*e**2/(2*R))+k*Z*e**2/(2*R)))
