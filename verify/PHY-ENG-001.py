# M11-1: uniform plank entering rough region, pulled by hanging mass via pulley.
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,M,g,L,mu,mu0,x = sp.symbols('m M g L mu mu0 x', positive=True)
v2 = (2*m*g*x - mu*M*g*x**2/L)/(M+m)        # energy: mgx - int_0^x mu M g x'/L dx' = (M+m)v^2/2
xs = sp.solve(sp.diff(v2,x),x)[0]; vmax = sp.sqrt(sp.simplify(v2.subs(x,xs)))
xmax = [r for r in sp.solve(v2,x) if r!=0][0]
print('(1) v^2=',sp.factor(v2))
print('(2) x at vmax',xs,' vmax=',sp.simplify(vmax),' xmax=',xmax,' cond xmax<=L -> mu>=2m/M')
# EOM approach (independent): (M+m) a = m g - mu M g x/L  -> SHM
om = sp.sqrt(mu*M*g/((M+m)*L)); Tstop = sp.pi/om
print('(3) t_stop=',sp.simplify(Tstop))
# (4) static friction needed at rest = m g ; available mu0 M g xmax/L
avail = sp.simplify(mu0*M*g*xmax/L); print('(4) available static friction', avail, '>= 2 m g >= m g since mu0>=mu')
# numeric check
vals={m:1.0,M:5.0,g:9.8,L:2.0,mu:0.6}
mv,Mv,gv,Lv,muv=[vals[s] for s in (m,M,g,L,mu)]
def f(t,y):
    xx,v=y; fr = muv*Mv*gv*min(max(xx,0),Lv)/Lv
    return [v,(mv*gv - fr*np.sign(v) if v!=0 else mv*gv-fr)/(Mv+mv)]
def stop(t,y): return y[1] if t>1e-6 else 1
stop.terminal=True; stop.direction=-1
r=solve_ivp(f,[0,50],[0,0],events=stop,rtol=1e-11,atol=1e-12,max_step=1e-3)
ts=r.t_events[0][0]; xe=r.y_events[0][0][0]
print('num t_stop',ts,' sym',float(Tstop.subs(vals)))
print('num xmax',xe,' sym',float(xmax.subs(vals)),' vmax num',r.y[1].max(),' sym',float(vmax.subs(vals)))
assert abs(ts-float(Tstop.subs(vals)))<1e-5 and abs(xe-float(xmax.subs(vals)))<1e-6
assert abs(r.y[1].max()-float(vmax.subs(vals)))<1e-4
print('ALL OK')
