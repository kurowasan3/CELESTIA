# M02-1: ball vs approaching smooth wall
import sympy as sp, math, numpy as np
from scipy.integrate import solve_ivp
v0,th,g,e,l,V,m=sp.symbols('v0 theta g e l V m',positive=True)
u=v0*sp.cos(th); w=v0*sp.sin(th)
t1=l/(u+V); T=2*w/g
ux_after=e*(u+V)+V   # speed away from wall (toward O)
Vs=sp.solve(sp.Eq(ux_after,u),V)[0]; print('V for equal speed:',sp.simplify(Vs))
# landing at O: x at collision u*t1, return u*(T-t1) -> t1=T/2
ls=sp.solve(sp.Eq((u*t1).subs(V,Vs),(ux_after*(T-t1)).subs(V,Vs)),l); print('l:',[sp.simplify(sp.trigsimp(x)) for x in ls])
J=m*(1+e)*(u+V); W=sp.simplify((J*V).subs(V,Vs)); Qh=sp.simplify((sp.Rational(1,2)*m*(1-e**2)*(u+V)**2).subs(V,Vs))
print('W',sp.factor(W),' Q',sp.factor(Qh),' diff',sp.simplify(W-Qh))
# Direct check: ground-frame KE before/after
KEb=sp.Rational(1,2)*m*(u**2); KEa=sp.Rational(1,2)*m*(ux_after.subs(V,Vs))**2; print('dKE horiz',sp.simplify(KEa-KEb))
# numeric event simulation
def sim(v0v,thv,gv,ev,lv,Vv):
    uu=v0v*math.cos(thv); ww=v0v*math.sin(thv)
    def f(t,y): return [y[2],y[3],0,-gv]
    def hitwall(t,y): return (lv-Vv*t)-y[0]
    hitwall.terminal=True; hitwall.direction=-1
    def floor(t,y): return y[1]
    floor.terminal=True; floor.direction=-1
    s1=solve_ivp(f,[0,10],[0,0,uu,ww],events=[hitwall,floor],rtol=1e-12,atol=1e-12)
    tc=s1.t[-1]; y=s1.y[:,-1]
    vx_rel=y[2]+Vv; y[2]=-ev*vx_rel-Vv
    s2=solve_ivp(f,[tc,10],y,events=[floor,hitwall],rtol=1e-12,atol=1e-12)
    return tc,s2.t[-1],s2.y[0,-1],y[2],uu
v0v,thv,gv,ev=10,0.7,9.8,0.6
uu=v0v*math.cos(thv); Vv=uu*(1-ev)/(1+ev); lv=v0v**2*math.sin(2*thv)/(gv*(1+ev))
print('numeric: tc, t_land, x_land, vx_after, u =',sim(v0v,thv,gv,ev,lv,Vv), ' T/2=',v0v*math.sin(thv)/gv)
