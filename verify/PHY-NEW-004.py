# M10-1: box Q (smooth floor) pushed by F=kt, block P inside with friction mu0/mu.
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,M,g,k,mu0,mu,t,tau = sp.symbols('m M g k mu0 mu t tau', positive=True)
t1 = sp.solve(sp.Eq(m*k*t/(m+M), mu0*m*g), t)[0]
print('(1) t1=',t1)
aP_after = mu*g
aQ_after = sp.simplify((k*t1 - mu*m*g)/M)
print('(2) before: both', sp.simplify(k*t1/(m+M)), ' after P', aP_after, ' Q', sp.factor(aQ_after))
alpha = (k*(t1+tau) - mu*m*g)/M - mu*g   # relative accel Q wrt P (Q - P)
u = sp.integrate(alpha,(tau,0,t1)); s = sp.integrate(sp.integrate(alpha,(tau,0,tau)),(tau,0,t1))
u_exp = (m+M)**2*mu0*(3*mu0-2*mu)*g**2/(2*M*k)
s_exp = (m+M)**3*mu0**2*(4*mu0-3*mu)*g**3/(6*M*k**2)
print('(3) u=',sp.factor(u), sp.simplify(u-u_exp)==0)
print('(4) s=',sp.factor(s), sp.simplify(s-s_exp)==0, ' heat=', sp.factor(mu*m*g*s))
assert sp.simplify(u-u_exp)==0 and sp.simplify(s-s_exp)==0
# numeric stick-slip simulation (event based)
vals = dict(m=1.0,M=2.0,g=9.8,k=3.0,mu0=0.5,mu=0.3)
mv,Mv,gv,kv,m0,m1 = [vals[x] for x in ['m','M','g','k','mu0','mu']]
def stick(t,y):  # y=[xP,vP,xQ,vQ] stuck together
    a = kv*t/(mv+Mv); return [y[1],a,y[3],a]
def ev(t,y): return mv*kv*t/(mv+Mv) - m0*mv*gv
ev.terminal=True
r1 = solve_ivp(stick,[0,100],[0,0,0,0],events=ev,rtol=1e-11,atol=1e-12)
t1n = r1.t_events[0][0]; y1 = r1.y_events[0][0]
def slip(t,y):
    aP = m1*gv; aQ=(kv*t - m1*mv*gv)/Mv; return [y[1],aP,y[3],aQ]
r2 = solve_ivp(slip,[t1n,2*t1n],y1,rtol=1e-11,atol=1e-12,dense_output=True)
yf = r2.y[:,-1]
# check relative velocity stays one-signed (no re-stick)
ts=np.linspace(t1n,2*t1n,400); Y=r2.sol(ts); assert np.all(Y[3]-Y[1]>=-1e-12)
sub={m:mv,M:Mv,g:gv,k:kv,mu0:m0,mu:m1}
print('t1 num',t1n,' sym',float(t1.subs(sub)))
print('u num',yf[3]-yf[1],' sym',float(u_exp.subs(sub)))
print('s num',(yf[2]-yf[0])-(y1[2]-y1[0]),' sym',float(s_exp.subs(sub)))
assert abs(yf[3]-yf[1]-float(u_exp.subs(sub)))<1e-6
assert abs((yf[2]-yf[0])-(y1[2]-y1[0])-float(s_exp.subs(sub)))<1e-6
print('ALL OK')
