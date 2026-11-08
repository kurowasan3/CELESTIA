# N45-2: rod enters field strip 0<=x<=w, rails closed by coil L, R=0
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
M,L,B,d,w,v0,x=sp.symbols('M L B d w v0 x',positive=True)
om=B*d/sp.sqrt(M*L)
# energy: 1/2 M v^2 + 1/2 L I^2 = 1/2 M v0^2 with I = B d x/L
I=B*d*x/L
v2=v0**2-L*I**2/M
print('v^2(x) =',sp.simplify(v2),' (om^2 x^2 form:',sp.simplify(v2-(v0**2-om**2*x**2)),')')
vth=sp.solve(sp.Eq(v2.subs(x,w),0),v0); print('threshold',vth)
v1=sp.sqrt(v2.subs(x,w)); print('v1',v1)
tin=sp.asin(om*w/v0)/om; print('time in field',tin)
# numeric
p=dict(M=0.3,L=0.05,B=0.6,d=0.4,w=0.5)
omn=p['B']*p['d']/np.sqrt(p['M']*p['L'])
def run(v0n,wall=None,T=20):
    # state x, v, I ; field strip [0,w]
    def f(t,y):
        X,V,Ii=y; inside=1.0 if 0<=X<=p['w'] else 0.0
        return [V,-inside*p['B']*p['d']*Ii/p['M'],inside*p['B']*p['d']*V/p['L']]
    evs=[]
    def exitw(t,y): return y[0]-p['w']
    def exit0(t,y): return y[0]+1e-9
    exit0.terminal=True; exit0.direction=-1
    s=solve_ivp(f,[0,T],[1e-12,v0n,0],events=[exitw,exit0],rtol=1e-11,atol=1e-13,max_step=1e-3)
    return s
vthn=p['B']*p['d']*p['w']/np.sqrt(p['M']*p['L']); print('vth num',vthn)
# pass-through case
v0n=1.5*vthn
s=run(v0n,T=3)
te=s.t_events[0][0]; ye=s.y_events[0][0]
print('exit time',te,'formula',np.arcsin(omn*p['w']/v0n)/omn)
print('exit v',ye[1],'formula',np.sqrt(v0n**2-omn**2*p['w']**2))
print('exit I',ye[2],'formula',p['B']*p['d']*p['w']/p['L'])
# rebound: start at x=w with -v1 and I=I1, run until x=0
def f2(t,y):
    X,V,Ii=y; inside=1.0 if 0<=X<=p['w'] else 0.0
    return [V,-inside*p['B']*p['d']*Ii/p['M'],inside*p['B']*p['d']*V/p['L']]
def at0(t,y): return y[0]
at0.terminal=True
s2=solve_ivp(f2,[0,5],[p['w']-1e-12,-ye[1],ye[2]],events=at0,rtol=1e-11,atol=1e-13,max_step=1e-3)
print('return to x=0: v',s2.y_events[0][0][1],'I',s2.y_events[0][0][2],' expect',-v0n,0)
# reflect case
v0n=0.6*vthn; s=run(v0n,T=5)
print('reflect: back at x=0 t',s.t_events[1][0],'formula pi/om',np.pi/omn,'v',s.y_events[1][0][1],'I',s.y_events[1][0][2])
