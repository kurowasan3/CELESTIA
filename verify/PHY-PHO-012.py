# J143-1 検証: 光電管の電流-電圧特性の形（放出方向の分布と電極の形）
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
rng=np.random.default_rng(1)
Nmc=400000
# (a) 立体角一様: 1-cos(theta) が [0,1] 一様 -> cos = 1-u
ca=1-rng.random(Nmc)
# (b) Lambert: sin^2 が一様 -> cos = sqrt(1-u)
cb=np.sqrt(1-rng.random(Nmc))
for x in [0.0,0.1,0.25,0.5,0.81,1.0]:   # x = eV/K0
    fa=np.mean(ca**2>=x); fb=np.mean(cb**2>=x)
    print(f'x={x}: MC(a)={fa:.4f} theory={1-np.sqrt(x):.4f} | MC(b)={fb:.4f} theory={1-x:.4f}')
# half current voltages
print('half (a): x=',0.25,' (b):',0.5)
# (5) 同心球: k=a/b
def frac_sphere(x,k):
    s=(1-x)/k**2
    if s>=1: return 1.0
    return 1-np.sqrt(1-s)
k=0.5
sa=np.sqrt(1-ca**2)
for x in [0.5,0.75,0.8,0.9,0.99]:
    mc=np.mean(1-x-sa**2*k**2>=0)
    print(f'sphere k={k} x={x}: MC={mc:.4f} theory={frac_sphere(x,k):.4f}')
# 軌道の数値積分で到達判定を確認（m=1, K0=1, a=1, b=2, 電位は同心球）
def reaches(theta,x,a=1.0,b=2.0):
    K0=1.0; m=1.0; v0=np.sqrt(2*K0/m)
    # U(r) = x*K0*(1/a-1/r)/(1/a-1/b)  (電子の位置エネルギー)
    c=x*K0/(1/a-1/b)
    def f(t,y):
        r=np.hypot(y[0],y[1]); F=-c/r**2  # 外向き正: 電子は内向きに押し戻される
        return [y[2],y[3],F*y[0]/r/m,F*y[1]/r/m]
    y0=[a,0,v0*np.cos(theta),v0*np.sin(theta)]
    ev_out=lambda t,y: np.hypot(y[0],y[1])-b; ev_out.terminal=True; ev_out.direction=1
    ev_in=lambda t,y: np.hypot(y[0],y[1])-a*0.999; ev_in.terminal=True; ev_in.direction=-1
    s=solve_ivp(f,[0,100],y0,events=[ev_out,ev_in],rtol=1e-10,atol=1e-12)
    return len(s.t_events[0])>0
x=0.8
thc=np.arcsin(np.sqrt((1-x)/k**2))
for th in [thc-0.02,thc+0.02, 0.3, 1.2]:
    print('theta',round(th,3),'reach',reaches(th,x),'pred',th<=thc)
X=sp.symbols('x'); print('(a) dI/dV at x=0 ->', sp.limit(sp.diff(1-sp.sqrt(X),X),X,0,'+'))
