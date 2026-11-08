import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
k,q,Q,F0,m,v0=sp.symbols('k q Q F0 m v0',positive=True)
rs=sp.sqrt(k*q*Q/F0); U=k*q*Q/rs
print('U_s=',sp.simplify(U))
v1=sp.sqrt(v0**2-2*U/m)
u=v1/2
rmin=k*q*Q/(U+sp.Rational(1,2)*m*v1**2-sp.Rational(1,2)*(2*m)*u**2)
print('rmin=',sp.simplify(rmin), ' target', sp.simplify(rmin-4*k*q*Q/(2*sp.sqrt(F0*k*q*Q)+m*v0**2)))
va,vb=sp.symbols('va vb',real=True)
sol=sp.solve([sp.Eq(m*va+m*vb,m*v1),sp.Eq(m*va**2/2+m*vb**2/2,m*v1**2/2+U)],[va,vb],dict=True)
for s_ in sol: print({kk:sp.simplify(vv) for kk,vv in s_.items()})
# numeric: k=q=Q=m=1, F0=1 -> rs=1, U=1; v0=2
K=1.0;F=1.0;M=1.0;V0=2.0
rs_=np.sqrt(K/F)
def phase1(t,y):
    xa,va_=y; r=0-xa; return [va_,-K/r**2/M]
ev=lambda t,y: (0-y[0])-rs_; ev.terminal=True; ev.direction=-1
s1=solve_ivp(phase1,[0,1e3],[-200,-np.sqrt(V0**2+2*K/200)],events=ev,rtol=1e-11,atol=1e-12)
# start far: energy matched so speed at infinity = V0; note va positive toward B
s1=solve_ivp(phase1,[0,1e3],[-200,np.sqrt(V0**2-2*K/200)],events=ev,rtol=1e-11,atol=1e-12)
xa,va1=s1.y_events[0][0]; print('released? v1 num',va1,' theory',np.sqrt(V0**2-2*1/1))
def phase2(t,y):
    xa,va_,xb,vb_=y; r=xb-xa; f=K/r**2; return [va_,-f/M,vb_,f/M]
s2=solve_ivp(phase2,[0,2e4],[xa,va1,0,0],rtol=1e-11,atol=1e-12,max_step=1)
r=s2.y[2]-s2.y[0]
print('rmin num',r.min(),' theory',4/(2+4))
print('final va,vb num',s2.y[1][-1],s2.y[3][-1],' theory',(np.sqrt(2)-np.sqrt(6))/2,(np.sqrt(2)+np.sqrt(6))/2)
