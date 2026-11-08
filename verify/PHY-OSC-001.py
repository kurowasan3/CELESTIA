# M25-1: elastic cord from O at height H above smooth floor; k=4mg/H, l0=15H/16
import sympy as sp, numpy as np
m,g,H,w,l=sp.symbols('m g H omega l',positive=True)
k=4*m*g/H; l0=sp.Rational(15,16)*H
# On floor: cord length l, angle a from vertical: cos a = H/l, radius r = sqrt(l^2-H^2)
T=k*(l-l0)
r=sp.sqrt(l**2-H**2)
# horizontal: T*r/l = m r w^2 -> T = m l w^2
lsol=sp.solve(sp.Eq(T,m*l*w**2),l)[0]; print('l(w)=',sp.simplify(lsol))
N=sp.simplify(m*g-T.subs(l,lsol)*H/lsol); print('N=',N)
assert sp.simplify(N-m*(g-H*w**2))==0
# circular motion on floor with r>0 requires l>H
wmin2=sp.solve(sp.Eq(lsol,H),w)[0]**2; print('w_min^2=',sp.simplify(wmin2))
wc2=sp.solve(sp.Eq(N,0),w)[0]**2; print('w_c^2=',wc2,' l at lift=',sp.simplify(lsol.subs(w,sp.sqrt(wc2))),' r=',sp.simplify(r.subs(l,lsol).subs(w,sp.sqrt(wc2))))
# after lift-off: conical pendulum. depth z below O: T cos a = mg, T sin a = m (l sin a) w^2 -> T = m l w^2, z = l cos a = g/w^2
z=sp.symbols('z',positive=True)
w2=2*g/H
lA=sp.solve(sp.Eq(k*(l-l0),m*l*w2),l)[0]; zA=g/w2
print('w^2=2g/H: l=',lA,' depth=',zA,' height above floor=',H-zA, 'r=',sp.sqrt(lA**2-zA**2))
# check with direct force balance numeric
Hv,gv,mv=1.0,9.8,1.0
kv=4*mv*gv/Hv; l0v=15/16*Hv; wv=np.sqrt(2*gv/Hv)
from scipy.optimize import fsolve
def F(x):
    L,a=x
    Tn=kv*(L-l0v)
    return [Tn*np.cos(a)-mv*gv, Tn*np.sin(a)-mv*L*np.sin(a)*wv**2]
L,a=fsolve(F,[1.8,1.2]); print('numeric l=',L,' height=',Hv-L*np.cos(a))
assert abs(L-15/8)<1e-9 and abs(Hv-L*np.cos(a)-0.5)<1e-9
# limit w^2 -> k/m = 4g/H: l -> oo, depth -> g/(4g/H)=H/4, height -> 3H/4
print('limit w^2=k/m:', sp.limit(lsol.subs(w,sp.sqrt(wc2)),w,0) if False else '', 'depth', g/(k/m),'height',sp.simplify(H-g/(k/m)))
# consistency: at rest sits on floor? hanging length l0+mg/k
print('hanging length',sp.simplify(l0+m*g/k),'> H')
# check N>0 and l>H inside window and T>0
for x in np.linspace(0.2501,0.9999,50):
    d={w:sp.sqrt(x*g/H)}
    assert sp.simplify(lsol.subs(d)/H)>1 and sp.simplify(N.subs(d)/(m*g))>0
print('OK')
