# T20233-1: two soap-film caps on tube ends; pressure max at hemisphere; stability of unequal caps
import sympy as sp, numpy as np
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
h,b,s,Rr=sp.symbols('h b sigma R',positive=True)
# virtual work for a bubble: dp*4piR^2 dR = 2 sigma d(4 pi R^2)
dp=sp.solve(sp.Eq(sp.Symbol('p')*4*sp.pi*Rr**2,2*s*8*sp.pi*Rr),sp.Symbol('p'))[0]; print('bubble dp=',dp)
Rc=(b**2+h**2)/(2*h); P=sp.simplify(4*s/Rc); V=sp.pi*h*(3*b**2+h**2)/6
print('P(h)=',P,' dP/dh zero at',sp.solve(sp.diff(P,h),h),' Pmax=',sp.simplify(P.subs(h,b)))
# check cap volume formula by integration
z=sp.symbols('z'); Vint=sp.integrate(sp.pi*(Rc**2-(z-(h-Rc))**2),(z,0,h)); print('V check',sp.simplify(Vint-V))
slope=sp.simplify(sp.diff(P,h)/sp.diff(V,h)); print('dP/dV =',sp.factor(slope))
# equal pressure distinct caps -> hA hB = b^2
hA,hB=sp.symbols('h_A h_B',positive=True)
print(sp.factor(sp.numer(sp.together(P.subs(h,hA)-P.subs(h,hB)))))
f=lambda hh: float(slope.subs({b:1,s:1,h:hh}))
print('(4) P(1/2),P(2)=',P.subs({b:1,s:1,h:sp.Rational(1,2)}),P.subs({b:1,s:1,h:2}),' sA+sB=',f(0.5)+f(2.0))
# ---- (4) general classification at fixed total volume ----
sl=lambda hh: 16*s*(b**2-hh**2)/(sp.pi*(b**2+hh**2)**3)
assert sp.simplify(sl(h)-slope)==0
tot=sp.factor(sp.simplify(sl(h)+sl(b**2/h))); print('(4i) k(hA)+k(b^2/hA) =',tot)
assert sp.simplify(tot-16*s*(b**2-h**2)**2/(sp.pi*b**4*(b**2+h**2)**2))==0
xx=sp.symbols('x',positive=True)
F=sp.pi*b**3/6*(xx*(3+xx**2)+(3*xx**2+1)/xx**3)
assert sp.simplify(F-(V.subs(h,xx*b)+V.subs(h,b/xx)))==0
Fp=sp.factor(sp.diff(F,xx)); print("(4iii) F'(x) =",Fp)
assert sp.simplify(sp.diff(F,xx)-sp.pi*b**3/2*(1+xx**2)*(1-1/xx**4))==0
print('(4ii) 2*V(hemisphere)=',sp.simplify(2*V.subs(h,b)),' F(1)=',sp.simplify(F.subs(xx,1)))
# dynamics: flow proportional to pressure difference, volume conserved; total volume Vt (units pi b^3)
Vf=lambda x: x*(3+x**2)/6
Pn=lambda hh: 8*hh/(1+hh*hh); dVdh=lambda hh: np.pi*(1+hh*hh)/2
def rhs(t,y):
    a,bb=y; q=Pn(bb)-Pn(a)
    return [q/dVdh(a), -q/dVdh(bb)]
for Vt in [1.2, 1.3, 1.4, 1.6, 2.5]:
    xs=brentq(lambda x: 2*Vf(x)-Vt,0.01,5.0)
    so=solve_ivp(rhs,[0,2000],[xs-1e-3,xs+1e-3],rtol=1e-10,atol=1e-12)
    a_,b_=so.y[:,-1]
    print('Vt=%.2f (crit %.4f) sym h=%.4f -> final hA=%.4f hB=%.4f hA*hB=%.4f dp=%.4f'%(Vt,4/3,xs,a_,b_,a_*b_,Pn(a_)))
    if Vt<4/3: assert abs(a_-xs)<1e-6 and abs(b_-xs)<1e-6
    else: assert abs(a_*b_-1)<1e-6 and a_<1<b_
# pressure decreases along asymmetric branch as total volume increases
prev=4.0
for Vt in [1.4,1.6,2.0,3.0]:
    x=brentq(lambda x: Vf(x)+Vf(1/x)-Vt,1e-3,1-1e-9); pv=Pn(x); assert pv<prev; prev=pv
    print('  branch Vt=%.1f x=%.4f dp=%.4f'%(Vt,x,pv))
print('OK')
