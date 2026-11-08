# Q154-1: べき乗則の磁力とばねのつり合いの消失 -> n=d/s, 浮上高さ, 微小振動の周期
import numpy as np, sympy as sp
from scipy.optimize import brentq
from scipy.integrate import solve_ivp
k=4.0; d=0.020; s=0.0050; m=5.0e-3; g=9.8
# symbolic tangency
K,n,x,u,kk=sp.symbols('K n x u k',positive=True)
sol=sp.solve([sp.Eq(kk*u,K/x**n),sp.Eq(kk,n*K/x**(n+1))],[K,n],dict=True)
print('tangency:',sol)
n_=d/s; K_=k*s*d**n_
print('n=',n_,'K=',K_)
# quasi-static approach simulation: D = distance from A's natural position to B. Track stable equilibrium u(D)
def eqs(D):
    f=lambda uu: k*uu-K_/(D-uu)**n_
    # stable root: smallest u in (0, D)
    us=np.linspace(1e-9,D*0.999,200001); fv=f(us)
    idx=np.where(np.sign(fv[:-1])!=np.sign(fv[1:]))[0]
    return [brentq(f,us[i],us[i+1]) for i in idx]
Ds=np.linspace(0.04,0.024,1601); last=None
for D in Ds:
    r=eqs(D)
    if not r: print('snap at D=',D,' last u=',last[0],' gap=',Dprev-last[0]); break
    last=r; Dprev=D
h=(K_/(m*g))**(1/n_); print('h=',h)
T=2*np.pi*np.sqrt(h/(n_*g)); print('T small-osc=',T)
# numerical oscillation in tube
def rhs(t,y): return [y[1], K_/(y[0])**n_/m - g]
amp=0.02*h
solv=solve_ivp(rhs,[0,1],[h-amp,0],max_step=1e-4,rtol=1e-10,atol=1e-13,dense_output=True)
ts=np.linspace(0,1,200001); xs=solv.sol(ts)[0]
# find successive minima
mins=[ts[i] for i in range(1,len(ts)-1) if xs[i]<xs[i-1] and xs[i]<=xs[i+1]]
print('numeric period=',np.diff(mins).mean())
