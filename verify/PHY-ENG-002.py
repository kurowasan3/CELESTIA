# M12-1: mass M at string midpoint, masses m at ends over smooth pegs (separation 2a)
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,M,g,a,y = sp.symbols('m M g a y', positive=True)
r = sp.sqrt(a**2+y**2)
U = -M*g*y + 2*m*g*(r-a)
y0 = sp.solve(sp.diff(U,y),y); print('(1) y0=',y0)
y0e = M*a/sp.sqrt(4*m**2-M**2)
y1 = [s for s in sp.solve(sp.Eq(M*y, 2*m*(r-a)),y) if s!=0]; print('(2) y1=',[sp.simplify(s) for s in y1])
y1e = 4*m*M*a/(4*m**2-M**2)
assert any(sp.simplify(s-y1e)==0 for s in y1)
# force-balance alternative for y0: 2 m g cos(th) = M g, cos=y/r
assert abs(float((2*m*y/r - M).subs(y,y0e).subs({m:1.3,M:1.1,a:0.7})))<1e-12
q = sp.symbols('q',positive=True)  # q = M/m
ratio = sp.simplify((y1e/y0e).subs(M,q*m)); print('ratio y1/y0=',ratio)
qs = sp.solve(sp.Eq(ratio,4),q); print('(3) M/m=',qs)
sub={M:sp.sqrt(3)*m}
Y0=sp.simplify(y0e.subs(sub)); Y1=sp.simplify(y1e.subs(sub)); print('y0,y1=',Y0,Y1)
# (4) speed at y0: KE = 1/2 M v^2 + 2*1/2 m (v y/r)^2
v=sp.symbols('v',positive=True)
KE = sp.Rational(1,2)*M*v**2 + m*(v*y/r)**2
vs = sp.solve(sp.Eq(KE, -U).subs(sub).subs(y,Y0), v)[0]
print('(4) v=',sp.nsimplify(sp.simplify(vs)), float(vs.subs({g:1,a:1})))
vexp = 2*sp.sqrt((2*sp.sqrt(3)-3)*g*a/3)
assert abs(float((vs-vexp).subs({g:9.8,a:0.7,m:1})))<1e-12
# (5) A rise = r(y1)-a ; hanging initially l/2 - a -> need l/2 - a - (r1-a) > 0
r1 = sp.simplify(r.subs(y,Y1)); print('(5) r1=',r1,' l_min=',2*r1)
# numeric: Lagrangian EOM for y
mv,Mv,gv,av = 1.0, np.sqrt(3), 9.8, 1.0
def f(t,s):
    yy,vv=s; R=np.sqrt(av**2+yy**2); c=yy/R
    # L = 1/2 (M + 2m c^2) ydot^2 - U ; d/dt: (M+2m c^2) ydd + 2m c dc/dy ydot^2 - ... 
    dc = av**2/R**3
    meff = Mv+2*mv*c**2
    ydd = (Mv*gv - 2*mv*gv*c - 2*mv*c*dc*vv**2)/meff
    return [vv,ydd]
ev=lambda t,s: s[1] if t>1e-6 else 1
ev.terminal=True; ev.direction=-1
sol=solve_ivp(f,[0,20],[0,0],events=ev,rtol=1e-11,atol=1e-12,dense_output=True)
y1n=sol.y_events[0][0][0]; print('num y1',y1n,' sym',float(Y1.subs(a,1)))
assert abs(y1n-4*np.sqrt(3))<1e-6
# speed when passing y0 = sqrt3
ev2=lambda t,s: s[0]-np.sqrt(3); ev2.terminal=True
s2=solve_ivp(f,[0,20],[0,0],events=ev2,rtol=1e-11,atol=1e-12)
print('num v at y0',s2.y_events[0][0][1],' sym',float(vexp.subs({g:9.8,a:1})))
assert abs(s2.y_events[0][0][1]-float(vexp.subs({g:9.8,a:1})))<1e-6
print('ALL OK')
# tension check: T = m(g + rddot), rddot = c*ydd + a^2/R^3 * ydot^2 (A,C rise by r-a)
ts=np.linspace(0,sol.t_events[0][0],2000); S=sol.sol(ts)
Tmin=1e9
for yy,vv in S.T:
    R=np.sqrt(av**2+yy**2); c=yy/R; ydd=f(0,[yy,vv])[1]
    Tmin=min(Tmin, mv*(gv + c*ydd + av**2/R**3*vv**2))
print('min tension / mg =',Tmin/(mv*gv)); assert Tmin>0
print('TENSION OK')
