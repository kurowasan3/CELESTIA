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
# (5) hA=1/2, hB=3/2
Vf=lambda x: x*(3+x**2)/6
Vt=Vf(0.5)+Vf(1.5); print('total V/(pi b^3)=',Vt,' sym each',Vt/2,' hemisphere',Vf(1.0))
x=brentq(lambda x: Vf(x)+Vf(1/x)-Vt,0.3,0.99); print('hA=',x,' hB=',1/x,' P=',8*x/(1+x*x),' stable sum slope=',f(x)+f(1/x))
print('initial pressures',8*0.5/1.25,8*1.5/(1+2.25))
# dynamics: flow rate proportional to pressure difference, volume conserved
Pn=lambda hh: 8*hh/(1+hh*hh); dVdh=lambda hh: np.pi*(1+hh*hh)/2
def rhs(t,y):
    a,bb=y; q=Pn(bb)-Pn(a)      # flow B->A
    return [q/dVdh(a), -q/dVdh(bb)]
sol=solve_ivp(rhs,[0,200],[0.5,1.5],rtol=1e-10,atol=1e-12)
print('dynamic final',sol.y[:,-1])
# perturb symmetric state with same total volume -> goes to asymmetric
xs=brentq(lambda x: 2*Vf(x)-Vt,1.0,2.0); sol2=solve_ivp(rhs,[0,400],[xs-1e-4,xs+1e-4],rtol=1e-10,atol=1e-12)
print('symmetric h=',xs,' perturbed ->',sol2.y[:,-1])
