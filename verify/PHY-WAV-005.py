# M68-1: wave reflected (fixed end) from a receding wall; reflected wave equation, beat period
import sympy as sp, numpy as np
A,f,v,u,L,t,x,tr=sp.symbols('A f v u L t x t_r',positive=True)
lam=v/f
y1=A*sp.sin(2*sp.pi*(f*t-x/lam))
X=L+u*t
fp=f*(v-u)/(v+u)
y2=-A*sp.sin(2*sp.pi*(fp*(t+(x-L)/v)-f*L/v))
# boundary condition at moving wall: y1+y2=0 at x=X(t)
bc=sp.simplify(sp.expand_trig((y1+y2).subs(x,X)))
print("BC residual:",sp.simplify(sp.trigsimp(bc)))
# check numerically too
fn=sp.lambdify((A,f,v,u,L,t),(y1+y2).subs(x,X))
print("BC num max:",max(abs(fn(1,3.0,10.0,1.3,4.0,tt)) for tt in np.linspace(0,20,2001)))
# y2 is function of t+x/v only (left-moving)
print("left-moving:",sp.simplify(sp.diff(y2,t)-v*sp.diff(y2,x)))
# left side: y3 = A sin 2pi f (t + x/v); total envelope
y3=A*sp.sin(2*sp.pi*f*(t+x/v))
# phase difference expressed with retarded wall position
tau=t+x/v
trr=(tau-L/v)/(1+u/v)
Xr=L+u*trr
ph3=2*sp.pi*f*tau
ph2=2*sp.pi*(fp*(t+(x-L)/v)-f*L/v)+sp.pi   # y2 = A sin(ph2)
print("phase diff - (4 pi Xr/lam + pi):",sp.simplify(ph3-ph2-(4*sp.pi*Xr/lam-sp.pi)))
print("beat freq:",sp.simplify(f-fp),"  dXr/dt:",sp.simplify(sp.diff(Xr,t)))
# numeric envelope period check
Av,fv,vv,uv,Lv=1.0,50.0,10.0,0.5,3.0
ts=np.linspace(0,10,400001); xo=-1.0
Y=(sp.lambdify((t,),(y2+y3).subs({A:Av,f:fv,v:vv,u:uv,L:Lv,x:xo}),'numpy'))(ts)
# envelope via local max in windows of one period
win=int((1/fv)/(ts[1]-ts[0]))
env=np.array([abs(Y[i:i+win]).max() for i in range(0,len(Y)-win,win)])
te=ts[0:len(Y)-win:win][:len(env)]
from scipy.signal import hilbert
envh=np.abs(hilbert(Y)); sl=slice(20000,-20000)
e=envh[sl]; tt=ts[sl]
idx=[i for i in range(1,len(e)-1) if e[i]<e[i-1] and e[i]<=e[i+1] and e[i]<0.3]
mins=[tt[i] for i in idx]
print("numeric envelope period:",np.diff(mins).mean()," predicted (v+u)/(2fu):",(vv+uv)/(2*fv*uv)," quasi-static lam/(2u):",vv/fv/(2*uv))
