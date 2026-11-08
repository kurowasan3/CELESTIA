import sympy as sp, numpy as np
t,v,g,k=sp.symbols('t v g k',positive=True)   # k = qE/(mg)
x=v*t-k*g*t**2/2; y=-g*t**2/2
r2=sp.expand(x**2+y**2); dr2=sp.expand(sp.diff(r2,t)/t)
print('dr2/t =',sp.collect(dr2,t))
disc=sp.discriminant(dr2,t); print('disc',sp.factor(disc))
print('roots k=3',sp.solve(dr2.subs(k,3),t))
# direction range for k=3: launch angle th
th=sp.symbols('theta')
ax,ay=-3*g,-g
def mono(thv):
    vv=1.0; gg=1.0
    vx,vy=np.cos(thv),np.sin(thv)
    tt=np.linspace(1e-4,20,20001)
    X=vx*tt-1.5*tt**2; Y=vy*tt-0.5*tt**2
    R=np.hypot(X,Y); return np.all(np.diff(R)>=-1e-12)
ths=np.radians(np.linspace(-10,50,60001))
ok=np.array([mono(a) for a in ths[::50]])
bad=np.degrees(ths[::50][~ok])
print('forbidden approx',bad.min(),bad.max())
print('pred',np.degrees(np.arctan(1/3)-np.arcsin(1/3)),np.degrees(np.arctan(1/3)+np.arcsin(1/3)))
