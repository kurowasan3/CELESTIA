# K20183-1 verification: isothermal gas column under piston / in closed box
import sympy as sp, numpy as np
N,m,g,k,T,M,S,z,H=sp.symbols('N m g k T M S z H',positive=True)
a=m*g/(k*T)
# (1) piston height
beta=N*m/M
h=sp.log(1+beta)/a
n0=N*a/(S*(1-sp.exp(-a*h)))
n=n0*sp.exp(-a*z)
assert sp.simplify(sp.integrate(n*S,(z,0,h))-N)==0
P=lambda zz:(n.subs(z,zz))*k*T
print('P(h)S-Mg =',sp.simplify(P(h)*S-M*g))
print('P(0)S-(M+Nm)g =',sp.simplify(P(0)*S-(M+N*m)*g))
# (2) gas PE
Ug=sp.integrate(m*g*z*n*S,(z,0,h))
print('Ug-(NkT-Mgh)=',sp.simplify(Ug-(N*k*T-M*g*h)))
# (3) heat capacity: E=3/2NkT+Ug+Mgh
E=sp.Rational(3,2)*N*k*T+Ug+M*g*h
print('dE/dT =',sp.simplify(sp.diff(E,T)))
print('d(Mgh)/dT =',sp.simplify(sp.diff(M*g*h,T)), ' expected Nk ln(1+b)/b')
# limits
b=sp.symbols('b',positive=True)
f=sp.log(1+b)/b
print('lim b->0',sp.limit(f,b,0),' lim b->oo',sp.limit(f,b,sp.oo))
# (5) closed box height H
nH=N*a/(S*(1-sp.exp(-a*H)))*sp.exp(-a*z)
UgH=sp.simplify(sp.integrate(m*g*z*nH*S,(z,0,H)))
x=sp.symbols('x',positive=True)
target=N*k*T*(1-x/(sp.exp(x)-1))
print('UgH check:',sp.simplify(UgH-target.subs(x,a*H)))
C=sp.diff(sp.Rational(3,2)*N*k*T+UgH,T)
Cform=N*k*(sp.Rational(5,2)-x**2*sp.exp(x)/(sp.exp(x)-1)**2)
print('C check:',sp.simplify(C-Cform.subs(x,a*H)))
cx=sp.Rational(5,2)-x**2*sp.exp(x)/(sp.exp(x)-1)**2
print('C/Nk x->0',sp.limit(cx,x,0),' x->oo',sp.limit(cx,x,sp.oo))
# numeric discrete-layer check of (3)
vals={N:1e3,m:1.0,g:1.0,k:1.0,M:500.0,S:1.0}
def Etot(TT):
    aa=1.0/TT; hh=np.log(1+2.0)/aa  # beta=1000/500=2
    zs=np.linspace(0,hh,200001); w=np.exp(-aa*zs); w/=np.trapezoid(w,zs)
    return 1.5*1e3*TT+1e3*np.trapezoid(zs*w,zs)+500*hh
d=1e-4; print('numeric C/Nk =',(Etot(1+d)-Etot(1-d))/(2*d)/1e3)
for xx in [0.5,2,5]:
    print('x',xx,'C/Nk',float(cx.subs(x,xx)))
# numeric checks of closed-box formulas (sympy simplify issue avoided)
import random
for _ in range(3):
    sub={N:random.uniform(1,5),m:random.uniform(.5,2),g:random.uniform(.5,2),k:1,T:random.uniform(.5,2),H:random.uniform(.2,4),S:1}
    xv=float((a*H).subs(sub))
    print('UgH diff',float((UgH-target.subs(x,a*H)).subs(sub)),' C diff',float((C-Cform.subs(x,a*H)).subs(sub)))
