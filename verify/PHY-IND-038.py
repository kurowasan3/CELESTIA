# O20232-1: コンデンサーの充放電電流がつくる磁場とコイルBの誘導電流（容量の2階微分）
import sympy as sp
t,a,b,d,er,e0,mu0,V,RA,RB,r,u,s0,w,K0=sp.symbols('t a b d epsilon_r epsilon_0 mu_0 V R_A R_B r u s_0 omega K0',positive=True)
Meff=mu0*sp.pi*RB**2/(2*RA)   # Φ_B = Meff * I_A
def IB_of_C(C):
    IA=V*sp.diff(C,t)               # 反時計回りを正とする充電電流
    emf=-Meff*sp.diff(IA,t)         # B の起電力（反時計回り正）
    return sp.simplify(IA), sp.simplify(emf/r)
# (1) 誘電体 s(t)
s=sp.Function('s')(t)
C1=e0*b*((a-s)+er*s)/d
IA1,IB1=IB_of_C(C1); print('(1) I_A=',IA1,'  I_B=',IB1)
print('    const speed:',sp.simplify(IB1.subs(s,u*t).doit()))
# (2) 振動
s2=a/2+s0*sp.sin(w*t)
IA2,IB2=IB_of_C(C1.subs(s,s2).doit()); print('(2) I_B=',sp.simplify(IB2))
amp=e0*b*(er-1)*V*Meff*w**2*s0/(d*r); print('    amp check:',sp.simplify(sp.Abs(IB2/sp.sin(w*t))-amp))
P=sp.simplify(r*amp**2/2); print('    avg power=',P)
# (3) 平行に一定の速さで離す
g=d+u*t; S=a*b
IA3,IB3=IB_of_C(e0*S/g); print('(3) I_A=',IA3,' I_B=',IB3)
Kc=sp.pi*RB**2*mu0*e0*S*V/(RA*r)
print('    check form:',sp.simplify(IB3 + Kc*u**2/g**3))
W=sp.integrate(r*IB3**2,(t,0,d/u)); print('    Joule heat=',sp.simplify(W), ' ratio to r Kc^2 u^3/d^5:',sp.simplify(W/(r*Kc**2*u**3/d**5)))
# (4) I_B=0 となる g(t)
G=sp.Function('G')(t)
ode=sp.diff(e0*S/G,t,2)
Y=sp.Function('Y')(t)  # Y=1/G
sol=sp.dsolve(sp.Eq(sp.diff(Y,t,2),0),Y,ics={Y.subs(t,0):1/d, sp.diff(Y,t).subs(t,0):-u/d**2})
print('(4) 1/g =',sol.rhs)
g4=d/(1-u*t/d)
print('    g(0)=',g4.subs(t,0),' gdot(0)=',sp.diff(g4,t).subs(t,0),' C2nd deriv=',sp.simplify(sp.diff(e0*S/g4,t,2)))
t2=sp.solve(sp.Eq(g4,2*d),t)[0]; print('    t at 2d=',t2,' speed=',sp.simplify(sp.diff(g4,t).subs(t,t2)))
# 数値検証：(3)のジュール熱を数値積分，(4)の I_B を差分で
import numpy as np
from scipy.integrate import quad
vals=dict(a=0.1,b=0.1,d=1e-3,u=1e-2,V=100.0,RA=0.5,RB=0.02,r=1.0,e0=8.85e-12,mu0=4e-7*np.pi)
Kn=np.pi*vals['RB']**2*vals['mu0']*vals['e0']*vals['a']*vals['b']*vals['V']/(vals['RA']*vals['r'])
def IBnum(tt,gfun,h=1e-6):
    C=lambda x: vals['e0']*vals['a']*vals['b']/gfun(x)
    d2C=(C(tt+h)-2*C(tt)+C(tt-h))/h**2
    return -np.pi*vals['RB']**2*vals['mu0']/(2*vals['RA'])*vals['V']*d2C/vals['r']
gl=lambda x: vals['d']+vals['u']*x
Wnum=quad(lambda x: vals['r']*IBnum(x,gl)**2,0,vals['d']/vals['u'])[0]
Wan=31/160*vals['r']*Kn**2*vals['u']**3/vals['d']**5
print('numeric W=',Wnum,' analytic=',Wan)
g4n=lambda x: vals['d']/(1-vals['u']*x/vals['d'])
print('numeric I_B for (4):',[IBnum(x,g4n,1e-5) for x in [0.01,0.03,0.05]], ' vs (3) scale',Kn*vals['u']**2/vals['d']**3)
