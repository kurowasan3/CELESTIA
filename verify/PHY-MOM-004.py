# M15-1: container against wall, spring+P inside; P released from compression a.
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,M,k,a = sp.symbols('m M k a', positive=True)
v1 = a*sp.sqrt(k/m)
u = m*v1/(m+M)
s = sp.sqrt(sp.simplify((m*v1**2 - (m+M)*u**2)/k)); print('(2) s=',sp.simplify(s))
assert sp.simplify(s - a*sp.sqrt(M/(m+M)))==0
V,vP = sp.symbols('V vP')
sol = sp.solve([sp.Eq(m*vP+M*V, m*v1), sp.Eq(m*vP**2+M*V**2, m*v1**2)],[vP,V],dict=True)
print('(3) natural-length states:',sol)
Vmax = 2*m*v1/(m+M); assert any(sp.simplify(d[V]-Vmax)==0 for d in sol)
print('(4) vP min = (m-M)/(m+M) v1 -> negative iff M>m')
# numeric with wall constraint
mv,Mv,kv,av = 1.0,3.0,10.0,0.2
def f(t,y):
    xP,vPv,X,Vv = y  # spring extension = xP - X
    F = -kv*(xP-X)    # force on P
    FQ = -F           # force on container from spring
    if X<=0 and Vv<=0 and FQ<0: FQ=0.0; aQ=0.0   # wall pushes
    else: aQ=FQ/Mv
    return [vPv,F/mv,Vv,aQ]
r=solve_ivp(f,[0,20],[-av,0,0,0],rtol=1e-10,atol=1e-12,max_step=1e-3)
xP,vPv,X,Vv=r.y
ext=xP-X
sub={m:mv,M:Mv,k:kv,a:av}
print('max ext num',ext.max(),' sym',float(s.subs(sub)))
print('Vmax num',Vv.max(),' sym',float(Vmax.subs(sub)))
print('vP min num',vPv[r.t>1].min(),' sym',float(((m-M)/(m+M)*v1).subs(sub)))
print('min container V after leaving',Vv.min(),' min X',X.min())
assert abs(ext.max()-float(s.subs(sub)))<1e-4 and abs(Vv.max()-float(Vmax.subs(sub)))<1e-4
assert abs(vPv[r.t>1].min()-float(((m-M)/(m+M)*v1).subs(sub)))<1e-4 and Vv.min()>-1e-6
print('ALL OK')
