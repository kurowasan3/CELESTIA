# S0101-1 verification: parachutist leaving horizontally-moving plane, linear (vector) drag after opening
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
g,h0,s,u,vT,w,V,X1,H = sp.symbols('g h0 s u v_T w V X1 H', positive=True)
k_over_m = g/vT
# exact relation (impulse of drag = -k * displacement), after opening:
t0 = sp.sqrt(2*s/g); vy0 = sp.sqrt(2*g*s)
# vertical (down +): m(V - vy0) = m g T - k (h0 - s)
T = sp.solve(sp.Eq(V - vy0, g*sp.Symbol('T') - k_over_m*(h0-s)), sp.Symbol('T'))[0]
# horizontal: m(w - u) = -k X_after
Xa = sp.solve(sp.Eq(w-u, -k_over_m*X1), X1)[0]
tL_exact = sp.simplify(t0 + T); X_exact = sp.simplify(u*t0 + Xa)
print('tL exact =', tL_exact); print('X exact =', X_exact)
tL_appr = sp.simplify(tL_exact.subs({V:vT})); X_appr = sp.simplify(X_exact.subs({w:0}))
print('tL approx =', tL_appr); print('X approx =', X_appr)
D = sp.symbols('D', positive=True)
s_sol = sp.solve(sp.Eq(X_appr, D), s); print('s for target D:', s_sol)
# numeric check with solve_ivp
gv=9.8; vTv=6.0; uv=40.0; h0v=1200.0; sv=300.0
km=gv/vTv
t0v=np.sqrt(2*sv/gv)
def f(t,y):
    x,yy,vx,vy=y
    return [vx,vy,-km*vx, gv-km*vy]
hit=lambda t,y: y[1]-h0v; hit.terminal=True; hit.direction=1
sol=solve_ivp(f,[t0v,500],[uv*t0v, sv, uv, gv*t0v],events=hit,rtol=1e-11,atol=1e-11)
tl=sol.t_events[0][0]; xl=sol.y_events[0][0][0]; vxl,vyl=sol.y_events[0][0][2:]
print('numeric tL, X, landing v =',tl,xl,vxl,vyl)
sub={g:gv,vT:vTv,u:uv,h0:h0v,s:sv}
print('exact formula tL (with measured V):',float(tL_exact.subs(sub).subs(V,vyl)),' X:',float(X_exact.subs(sub).subs(w,vxl)))
print('approx tL:',float(tL_appr.subs(sub)),' approx X:',float(X_appr.subs(sub)))
# target D check
Dv=600.0; sD=float(s_sol[0].subs({g:gv,vT:vTv,u:uv,D:Dv})); print('s for D=600:',sD, ' formula', gv/2*(Dv/uv-vTv/gv)**2)
t0v=np.sqrt(2*sD/gv)
sol=solve_ivp(f,[t0v,500],[uv*t0v, sD, uv, gv*t0v],events=hit,rtol=1e-11,atol=1e-11)
print('numeric landing x for that s:',sol.y_events[0][0][0])
