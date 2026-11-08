# O20091-1 : 床と衝突するばね振り子の衝突間隔から e と床の位置を逆算
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
w,d,h,u,e=sp.symbols('omega d h u e',positive=True)
T=2*sp.pi/w
dt=lambda uu: 2/w*(sp.pi-sp.atan(uu/(w*d)))
u1=sp.solve(sp.Eq(dt(u),2*T/3),u); u2=sp.solve(sp.Eq(dt(u),3*T/4),u); print('u1',u1,' u2',u2)
ev=sp.simplify(u2[0]/u1[0]); print('e=',ev)
u0=u1[0]/ev; hs=sp.solve(sp.Eq(w**2*h**2,u0**2+w**2*d**2),d); print('d=',hs)
print('third interval /T:',sp.simplify(dt(u2[0]*ev)/T))
t1=(sp.pi-sp.acos(d/h))/w; print('t1 with d=h/sqrt10:',sp.simplify(t1.subs(d,h/sp.sqrt(10))), float((sp.pi-sp.atan(3))))
# 数値シミュレーション
wv,hv=1.0,1.0; dv=hv/np.sqrt(10); evv=1/np.sqrt(3)
def f(t,y): return [y[1],-wv**2*y[0]]
def hit(t,y): return y[0]+dv
hit.terminal=True; hit.direction=-1
y=[hv,0]; t0=0; times=[]
for i in range(8):
    s=solve_ivp(f,(t0,t0+50),y,events=hit,rtol=1e-12,atol=1e-13)
    tc=s.t_events[0][0]; yc=s.y_events[0][0]; times.append(tc)
    y=[-dv+1e-12, -evv*yc[1]]; t0=tc
Tn=2*np.pi/wv
print('t1=',times[0],' theory',(np.pi-np.arctan(3))/wv)
print('intervals/T:',[round((b-a)/Tn,6) for a,b in zip(times,times[1:])])
print('energy lost frac (h^2-d^2)/h^2 =',(hv**2-dv**2)/hv**2)
