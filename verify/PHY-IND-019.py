# T20092-1 verification: coil through two adjacent opposite-field strips
import numpy as np
from scipy.integrate import solve_ivp
m,R,B,b,a,g=0.05,0.1,2.0,0.2,0.02,9.8
K=B**2*b**2/R
v0=m*g/K; h=v0**2/(2*g)
def fac(x):  # x: displacement of lower side below P
    if 0<=x<a: return 1
    if a<=x<2*a: return 4
    if 2*a<=x<3*a: return 1
    return 0
def rhs(t,y):
    x,v=y; return [v, g-fac(x)*K*v/m]
def ev(X):
    f=lambda t,y: y[0]-X; f.terminal=True; f.direction=1; return f
# start at P with speed v0 (released from height h)
# check free fall from h gives v0
print('v0',v0,'sqrt(2gh)',np.sqrt(2*g*h))
res={}
y=[0,v0]; t=0; heat=0
for X in [a,2*a,3*a]:
    s=solve_ivp(rhs,[t,t+50],y,events=ev(X),max_step=1e-4,rtol=1e-10,atol=1e-12,dense_output=True)
    tt=s.t_events[0][0]; yy=s.y_events[0][0]
    # heat = integral of fac*K*v^2
    ts=np.linspace(t,tt,20001); vs=s.sol(ts)[1]; xs=s.sol(ts)[0]
    fs=np.array([fac(xx-1e-12) if xx>=X-1e-9 else fac(xx) for xx in xs])
    q=np.trapezoid(fs*K*vs**2,ts)
    res[X]=(tt-t,yy[1],q); t=tt; y=list(yy)
t1,vA,q1=res[a]; t2,u,q2=res[2*a]; t3,w,q3=res[3*a]
print('phase1 time',t1,'a/v0',a/v0,'v end',vA)
print('u',u,'> v0/4',v0/4)
print('t2',t2,'formula',4*a/v0-(v0-u)/g)
print('Q2',q2,'formula',m*g*a+0.5*m*(v0**2-u**2))
print('t3',t3,'formula',a/v0+(w-u)/g, 'w<v0',w,v0)
print('total',t1+t2+t3,'formula',6*a/v0+(w-v0)/g)
