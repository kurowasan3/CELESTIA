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
# continuity check of energy at r=R between models? E_out(r)= -kZe^2/(2r); E_in(r)= m w^2 r^2 - 3kZe^2/(2R); at r=R both -kZe^2/(2R)
print(sp.simplify((m*w**2*R**2-3*k*Z*e**2/(2*R))+k*Z*e**2/(2*R)))

# ---- new (4): self-consistency equivalence and energy shift ----
rho,RR,y,K=sp.symbols('rho R y K',positive=True)
# s^4 = rho R^3  (from r_n^4 = n^2 a0 R^3/(mu Z) and rho = n^2 a0/(mu Z))
print('s^4 - rho R^3 =',sp.simplify(rn4 - (n**2*a0/(mu*Z))*R**3))
for rv in [0.3,0.9,1.0,1.7,5.0]:
    sv=(rv*1.0**3)**0.25   # R=1
    assert (sv<1.0)==(rv<1.0) or abs(rv-1)<1e-12
# energies with K=kZe^2, m w^2 = K/R^3
s_=(rho*RR**3)**sp.Rational(1,4)
E3=K/RR**3*s_**2-3*K/(2*RR); E1=-K/(2*rho)
delta=sp.simplify((E3-E1).subs(rho,RR/y**2))
target=K/RR*(y-1)**2*(y+2)/(2*y)
print('delta - target =',sp.simplify(delta-target))
# K = 2 Z EH a0 -> delta = Z EH a0/R (y-1)^2(y+2)/y
print('K check (k e^2=2EH a0):',sp.simplify((k*Z*e**2).subs(subs)-2*Z*EH*a0))
# also check E3 directly from sympy model expression En (n hbar w - 3kZe^2/2R) with rho=n^2 hbar^2/(m k Z e^2)
rho_expr=n**2*hb**2/(m*k*Z*e**2)
d_direct=sp.simplify(En - (-k*Z*e**2/(2*rho_expr)))
yv=sp.sqrt(R/rho_expr)
print('direct delta - formula =',sp.simplify(d_direct - (k*Z*e**2/R*(yv-1)**2*(yv+2)/(2*yv))))
f=sp.lambdify(y,(y-1)**2*(y+2)/y)
ys=np.linspace(0.01,10,2000); assert (f(ys)>=-1e-15).all(); print('delta>=0 on grid; delta(1)=',f(1.0))
print('series near y=1:',sp.series((y-1)**2*(y+2)/y,y,1,3))
