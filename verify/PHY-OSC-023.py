# T20191-1: inverted pendulum on cart, transport by three half-period SHM segments
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
t,g,l,th1=sp.symbols('t g ell theta1',positive=True)
w=sp.sqrt(g/l); Th=sp.pi/w
# segments (local time s in [0,Th])
seg=[th1/2*(1-sp.cos(w*t)), th1*sp.cos(w*t), -th1/2*(1+sp.cos(w*t))]
# check continuity of theta and theta'
for i in range(2):
    print('cont',sp.simplify(seg[i].subs(t,Th)-seg[i+1].subs(t,0)), sp.simplify(sp.diff(seg[i],t).subs(t,Th)-sp.diff(seg[i+1],t).subs(t,0)))
print('start/end', seg[0].subs(t,0), sp.diff(seg[0],t).subs(t,0), sp.simplify(seg[2].subs(t,Th)), sp.simplify(sp.diff(seg[2],t).subs(t,Th)))
# cart acceleration a = g th - l th''
acc=[sp.simplify(g*s-l*sp.diff(s,t,2)) for s in seg]; print('a(t) segments:',acc)
print('a(0+)=',acc[0].subs(t,0))
# cart velocity/position by integration
v=0; x=0; X=[];V=[]
for a in acc:
    vt=v+sp.integrate(a,(t,0,t)); xt=x+sp.integrate(vt,(t,0,t))
    V.append(sp.simplify(vt)); X.append(sp.simplify(xt))
    v=sp.simplify(vt.subs(t,Th)); x=sp.simplify(xt.subs(t,Th))
print('final v=',v,' final x=',sp.simplify(x))
print('v at end seg1:',sp.simplify(V[0].subs(t,Th)))
# max speed in seg2
vs2=V[1]; print('v seg2:',vs2, ' at Th/2:',sp.simplify(vs2.subs(t,Th/2)))
print('integral theta total:',sp.simplify(sum(sp.integrate(s,(t,0,Th)) for s in seg)))
# numeric: simulate nonlinear-free linear system with prescribed a(t) and check theta follows, and cart speed max
gv,lv,t1=9.8,1.0,0.05; wv=np.sqrt(gv/lv); Tv=np.pi/wv
def afun(tt):
    k=int(min(tt//Tv,2)); s=tt-k*Tv
    return [gv*t1/2*(1-2*np.cos(wv*s)), 2*gv*t1*np.cos(wv*s), -gv*t1/2*(1+2*np.cos(wv*s))][k]
# check acc expressions numerically equal sympy
f=[sp.lambdify(t,a.subs({g:gv,l:lv,th1:t1})) for a in acc]
print([ (f[k](0.3), afun(k*Tv+0.3)) for k in range(3)])
def rhs(tt,y):
    th,om,xc,vc=y; a=afun(tt)
    return [om,(gv*th-a)/lv,vc,a]
sol=solve_ivp(rhs,[0,3*Tv],[0,0,0,0],max_step=1e-3,rtol=1e-10,atol=1e-12)
print('numeric final theta,omega,x,v:',sol.y[:,-1], ' pred x=',np.pi**2*lv*t1)
print('numeric max v:',sol.y[3].max(),' pred:',gv*t1/wv*(np.pi/2+2), ' min v:',sol.y[3].min())
# nonlinear check with exact inverted pendulum eq l th''= g sin th - a cos th
def rhsn(tt,y):
    th,om,xc,vc=y; a=afun(tt)
    return [om,(gv*np.sin(th)-a*np.cos(th))/lv,vc,a]
soln=solve_ivp(rhsn,[0,3*Tv],[0,0,0,0],max_step=1e-3,rtol=1e-10,atol=1e-12)
print('nonlinear final:',soln.y[:,-1])
