# T20151-1: release of the pivot of a pendulum -> CM free fall + uniform rotation
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,M,l,g,v0,al=sp.symbols('m M l g v0 alpha',positive=True)
# (1) slack angle with B fixed
v2=v0**2-2*g*l*(1+sp.sin(al))
T=m*v2/l-m*g*sp.sin(al)
print('(1) sin alpha =',sp.solve(sp.Eq(T,0),sp.sin(al)))
u2=v0**2-2*g*l
T1=m*u2/l; mu=m*M/(m+M); T2=mu*u2/l
print('(2) T1,T2,ratio',T1,sp.simplify(T2),sp.simplify(T2/T1))
u=sp.sqrt(u2); t1=sp.pi*l/(2*u)
yB=m*u*t1/(m+M)-g*t1**2/2-m*l/(m+M)
print('(5) yB=',sp.simplify(sp.expand(yB)))
# numeric check with stiff spring
vals=dict(m=1.0,M=2.0,l=1.0,g=9.8); vv0=np.sqrt(3.5*9.8*1.0)
mm,MM,ll,gg=vals['m'],vals['M'],vals['l'],vals['g']
uu=np.sqrt(vv0**2-2*gg*ll)
K=2e7
def f(t,s):
    xa,ya,xb,yb,vxa,vya,vxb,vyb=s
    dx,dy=xa-xb,ya-yb; r=np.hypot(dx,dy); Ten=K*(r-ll)
    ex,ey=dx/r,dy/r
    return [vxa,vya,vxb,vyb,-Ten*ex/mm,-gg-Ten*ey/mm,Ten*ex/MM,-gg+Ten*ey/MM]
# static stretch for centripetal at t=0 to reduce ringing
Ttheory=mm*MM/(mm+MM)*uu**2/ll
s0=[ll+Ttheory/K,0,0,0,0,uu,0,0]
tt1=np.pi*ll/(2*uu)
sol=solve_ivp(f,[0,tt1],s0,rtol=1e-10,atol=1e-12,method='DOP853',dense_output=True)
xa,ya,xb,yb=sol.y[:4,-1]
print('num A rel B at t1 (should be ~(0,l)):',xa-xb,ya-yb)
yBth=float(yB.subs({m:mm,M:MM,l:ll,g:gg,v0:vv0}))
print('num yB',yb,' theory',yBth,' xB num',xb,' theory',mm*ll/(mm+MM))
ts=np.linspace(0,tt1,50); S=sol.sol(ts)
Tn=K*(np.hypot(S[0]-S[2],S[1]-S[3])-ll)
print('tension num mean/min/max',Tn.mean(),Tn.min(),Tn.max(),' theory',Ttheory)
