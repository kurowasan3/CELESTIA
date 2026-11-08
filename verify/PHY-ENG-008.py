import sympy as sp
m,M,u,w=sp.symbols('m M u w',positive=True)
def state(vcar):
    Vp=-m*vcar/M   # platform velocity (ice frame)
    return vcar,Vp
def K(vc,vp,shift=0):
    return m*(vc+shift)**2/2+M*(vp+shift)**2/2
c0,p0=state(0); c1,p1=state(u); c2,p2=state(2*u)
E1=sp.simplify(K(c1,p1)-K(c0,p0)); E2=sp.simplify(K(c2,p2)-K(c1,p1))
print('E1',sp.factor(E1),'E2',sp.factor(E2),'ratio',sp.simplify(E2/E1))
# observer moving at -w (velocities shift by +w)
dKc=sp.expand(m*(c2+w)**2/2-m*(c1+w)**2/2); dKp=sp.expand(M*(p2+w)**2/2-M*(p1+w)**2/2)
print('dKcar',dKc,' dKplat',dKp,' sum-E2',sp.simplify(dKc+dKp-E2))
print('limit dKplat M->oo',sp.limit(dKp,M,sp.oo),' ice frame',sp.limit(M*p2**2/2-M*p1**2/2,M,sp.oo))
