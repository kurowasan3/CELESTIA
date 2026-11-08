import sympy as sp, numpy as np
m,g,R,v,L,th=sp.symbols('m g R v L theta',positive=True)
r=sp.Matrix([R*sp.cos(th),R*sp.sin(th)])
vel=sp.diff(r,th)*v/R          # counterclockwise, speed v
acc=sp.diff(vel,th)*v/R
T=sp.simplify(m*acc-m*sp.Matrix([0,-g]))
print('T vec',T.T)
print('|T|^2',sp.simplify(T.dot(T)))
P=sp.simplify(T.dot(vel)); print('power',P)
for t0 in [sp.pi/2,3*sp.pi/2,0]:
    Tv=T.subs(th,t0); n=sp.sqrt(Tv.dot(Tv))
    print(t0,'T=',sp.simplify(Tv.T),'hand=',sp.simplify((r.subs(th,t0)+L*Tv/n).T))
# numeric check: v^2<gR at top => hand above
vals={m:1,g:9.8,R:0.5,v:1.5,L:0.5}
Tv=T.subs(th,sp.pi/2).subs(vals); print('slow case top T',Tv.T)
print('cycle work',sp.integrate(P/v*R,(th,0,2*sp.pi)))
