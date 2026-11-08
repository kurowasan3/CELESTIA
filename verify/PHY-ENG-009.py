# Q066-1: paths from A to B with speed-preserving deflections
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
L,H,g=sp.symbols('L H g',positive=True)
T0=sp.sqrt(2*(L**2+H**2)/(g*H))       # straight smooth incline from rest
T1=sp.sqrt(2*H/g)+L/sp.sqrt(2*g*H)     # drop H then horizontal
print('T1 simplified',sp.simplify(T1-(2*H+L)/sp.sqrt(2*g*H)))
cond=sp.solve(sp.simplify(T1**2-T0**2)<0,L); print('T1<T0 when',cond)
num={L:2,H:1,g:sp.Rational(98,10)}
print('T0',sp.N(T0.subs(num)),'T1',sp.N(T1.subs(num)),'bound L/sqrt(2gH)',sp.N((L/sp.sqrt(2*g*H)).subs(num)))
# numerical check of straight incline time
th=np.arctan2(1,2); gn=9.8
r=solve_ivp(lambda t,y:[y[1],gn*np.sin(th)],[0,5],[0,0],events=lambda t,y:y[0]-np.sqrt(5),rtol=1e-11)
print('incline numeric',r.t_events[0])
# drop then horizontal numeric
v=np.sqrt(2*gn*1); print('L-path numeric',np.sqrt(2/gn)+2/v)
# random paths: piecewise straight smooth segments, check all times exceed bound
rng=np.random.default_rng(0); best=1e9
def ptime(ys):  # ys: depths at equally spaced x nodes, endpoints 0 and 1 ; segments straight, speed sqrt(2 g y)
    xs=np.linspace(0,2,len(ys)); T=0
    for i in range(len(ys)-1):
        ds=np.hypot(xs[i+1]-xs[i],ys[i+1]-ys[i]); v1=np.sqrt(2*gn*ys[i]); v2=np.sqrt(2*gn*ys[i+1])
        if v1+v2==0: return np.inf
        T+=2*ds/(v1+v2)   # uniform acceleration along straight segment
    return T
for _ in range(20000):
    n=rng.integers(2,8); ys=np.concatenate([[0],rng.uniform(0,1,n),[1]])  # board: depth <= H
    best=min(best,ptime(ys))
print('best random piecewise path (depth<=H)',best)
# vertical drop node first then optimized: scipy minimize
from scipy.optimize import minimize
def f(z):
    ys=np.concatenate([[0],np.clip(np.abs(z),0,1),[1]]); xs=np.concatenate([[0],np.linspace(0,2,len(z)+2)[1:-1],[2]])
    return ptime(ys)
res=minimize(f,np.full(30,1.0),method='Powell',options={'maxiter':200000,'xtol':1e-8,'ftol':1e-12}); print('optimized piecewise (depth<=H)',res.fun,' > bound?',res.fun>0.4518)
# cycloid (brachistochrone) exact time for reference
from scipy.optimize import brentq
phi=brentq(lambda p:(p-np.sin(p))/(1-np.cos(p))-2,0.1,2*np.pi-0.01); R=1/(1-np.cos(phi)); print('cycloid time',phi*np.sqrt(R/gn))
# staircase path: drop h, horizontal L, drop H-h  -> T(h)=sqrt(2H/g)+L/sqrt(2gh)
h=sp.symbols('h',positive=True)
Th=sp.sqrt(2*h/g)+L/sp.sqrt(2*g*h)+(sp.sqrt(2*g*H)-sp.sqrt(2*g*h))/g
print('T(h) - formula',sp.simplify(Th-(sp.sqrt(2*H/g)+L/sp.sqrt(2*g*h))))
def stair(hn):
    t1=np.sqrt(2*hn/gn); v=np.sqrt(2*gn*hn); t2=2/v
    r=solve_ivp(lambda t,y:[y[1],gn],[0,5],[0,0],events=lambda t,y:y[0]-(1-hn),rtol=1e-12)  # vertical drop from rest? no: start with speed v downward
    r=solve_ivp(lambda t,y:[y[1],gn],[0,5],[0,v],events=lambda t,y:y[0]-(1-hn),rtol=1e-12)
    return t1+t2+r.t_events[0][0]
for hn in [0.25,0.5,0.75]:
    print(hn,stair(hn),float((sp.sqrt(2*H/g)+L/sp.sqrt(2*g*h)).subs(num).subs(h,hn)))
