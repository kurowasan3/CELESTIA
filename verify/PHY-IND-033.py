# J127-1 検証：一定速度で駆動される棒が自由な棒を引きずる過程のジュール熱（段階的・連続的な加速）
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp, quad
t,m,B,l,R,w,al,T=sp.symbols('t m B l R w alpha T',positive=True)
tau=2*m*R/(B**2*l**2)
v=sp.Function('v')
sol=sp.dsolve(sp.Eq(m*v(t).diff(t),B*l*B*l*(w-v(t))/(2*R)),v(t),ics={v(0):0}).rhs
print('(1) v(t)=',sp.simplify(sol))
I=B*l*(w-sol)/(2*R)
H=sp.integrate(2*R*I**2,(t,0,sp.oo)); print('(2) heat=',sp.simplify(H))
Wext=sp.integrate(B*l*I*w,(t,0,sp.oo)); print('    Wext=',sp.simplify(Wext))
dx=sp.integrate(w-sol,(t,0,sp.oo)); print('(3) lag=',sp.simplify(dx))
# (4) N段階
N=sp.symbols('N',positive=True,integer=True)
print('(4) N steps heat =',sp.simplify(N*sp.Rational(1,2)*m*(w/N)**2))
# (5) 一定加速度 alpha で 0->w（時間 T=w/alpha）
s=sp.Function('s')
ss=sp.dsolve(sp.Eq(s(t).diff(t),al-s(t)/tau),s(t),ics={s(0):0}).rhs
k=B**2*l**2/(2*R)
H5=sp.integrate(k*ss**2,(t,0,w/al))+k*(ss.subs(t,w/al))**2*tau/2
H5=sp.simplify(H5.subs(al,w/T))
target=m*w**2*(tau/T)*(1-(tau/T)*(1-sp.exp(-T/tau)))
print('(5) diff=',sp.simplify(H5-target))
print('   T->0:',sp.limit(target,T,0),'  T->oo leading:',sp.limit(target*T,T,sp.oo))
# 数値積分で (5) 照合
mn,Bn,ln,Rn,wn,Tn=0.1,0.5,0.4,0.2,2.0,1.5
taun=2*mn*Rn/(Bn**2*ln**2); kn=Bn**2*ln**2/(2*Rn)
uq=lambda tt: wn*min(tt/Tn,1.0)
def f(tt,y): 
    rel=uq(tt)-y[0]; return [kn*rel/mn, kn*rel**2]
r=solve_ivp(f,[0,Tn+40*taun],[0,0],rtol=1e-10,atol=1e-12,max_step=0.01)
print('num heat=',r.y[1,-1],' formula=',float(target.subs({m:mn,w:wn,T:Tn,B:Bn,l:ln,R:Rn})))
# (2) 数値
r2=solve_ivp(lambda tt,y:[kn*(wn-y[0])/mn,kn*(wn-y[0])**2,(wn-y[0])],[0,60*taun],[0,0,0],rtol=1e-10,atol=1e-12)
print('num heat2=',r2.y[1,-1],' 0.5mw^2=',0.5*mn*wn**2,' lag=',r2.y[2,-1],' formula',float(dx.subs({m:mn,w:wn,B:Bn,l:ln,R:Rn})))
