# K20132-1 verification
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
Q,q,e0,S,d,dd,m,g,k,t=sp.symbols('Q q epsilon0 S d Delta m g k t',positive=True)
# (1) force at fixed charge via energy
W=lambda x: Q**2*x/(2*e0*S)
F=sp.simplify((W(d+dd)-W(d))/dd); print('(1) F =',F)
F1=Q**2/(2*e0*S); F2=q**2/(2*e0*S)
T0=sp.symbols('T0')
sol=sp.solve([sp.Eq(F1,m*g+T0), sp.Eq(T0,m*g+F2)],[Q,T0],dict=True)[0]
print('(2) Q=',sp.simplify(sol[Q]),' T0=',sol[T0],' ext=',sol[T0]/k)
# (3)(4) ODE: xA,xB displacements (up +) from initial, after B neutralized
w=sp.sqrt(2*k/m)
xA=sp.Function('xA'); xB=sp.Function('xB')
T=lambda: sol[T0]-k*(xA(t)-xB(t))   # tension (spring between, A above B): A up -> spring longer
# A above: spring extension increases when xB decreases or xA increases? spring hangs from A down to B: length = (yA - yB); extension up when xA - xB increases
Ten = sol[T0] + k*(xA(t)-xB(t))
eqs=[sp.Eq(m*xA(t).diff(t,2), F1.subs(Q,sol[Q]) - m*g - Ten), sp.Eq(m*xB(t).diff(t,2), Ten - m*g)]
s=sp.dsolve(eqs,[xA(t),xB(t)],ics={xA(0):0,xB(0):0,xA(t).diff(t).subs(t,0):0,xB(t).diff(t).subs(t,0):0})
xa=sp.simplify(s[0].rhs); xb=sp.simplify(s[1].rhs)
print('xA =',xa); print('xB =',xb)
claimA = F2*t**2/(4*m) - F2/(4*k)*(1-sp.cos(w*t))
print('check xA:',sp.simplify(sp.expand_complex(xa-claimA).rewrite(sp.cos)))
aA=sp.simplify(xa.diff(t,2)); print('aA=',sp.simplify(aA.rewrite(sp.cos)), ' claim (F2/2m)(1-cos):', sp.simplify((aA-(F2/(2*m))*(1-sp.cos(w*t))).rewrite(sp.cos)))
Tt=sp.simplify((sol[T0]+k*(xa-xb)).rewrite(sp.cos)); print('Tension(t)=',Tt)
Tp=2*sp.pi/w
print('(5) d=',sp.simplify(xa.subs(t,Tp)), ' vA=',sp.simplify(xa.diff(t).subs(t,Tp)))
# numeric check with solve_ivp
vals={m:0.2,g:9.8,k:50.0,e0:8.85e-12,S:0.01,q:3e-6}
F2n=float(F2.subs(vals)); T0n=float(sol[T0].subs(vals)); F1n=float(F1.subs(Q,sol[Q]).subs(vals))
mn,kn,gn=0.2,50.0,9.8
def f(tt,y):
    xa_,xb_,va,vb=y; Te=T0n+kn*(xa_-xb_)
    return [va,vb,(F1n-mn*gn-Te)/mn,(Te-mn*gn)/mn]
wn=np.sqrt(2*kn/mn); Tpn=2*np.pi/wn
r=solve_ivp(f,[0,Tpn],[0,0,0,0],rtol=1e-11,atol=1e-14,dense_output=True)
print('numeric xA(T)=',r.y[0,-1],' formula=',np.pi**2*F2n/(2*kn),' vA=',r.y[2,-1],' formula=',np.pi*F2n/np.sqrt(2*mn*kn))
ts=np.linspace(0,Tpn,2000); Y=r.sol(ts)
print('min vA over (0,T]=',Y[2,1:].min(),' min tension=',(T0n+kn*(Y[0]-Y[1])).min(),' mg=',mn*gn)
