# M26-1: bead on rod at angle th from vertical axis, rotating with w; static friction mu (0<mu<1); bead at distance s from O along rod.
import sympy as sp, numpy as np
m,g,s,w,mu,th=sp.symbols('m g s omega mu theta',positive=True)
r=s*sp.sin(th)
Fpar=m*r*w**2*sp.sin(th)-m*g*sp.cos(th)     # along rod, outward/up positive (rotating frame)
Nn=m*r*w**2*sp.cos(th)+m*g*sp.sin(th)        # magnitude of normal force needed (>0 always for 0<th<90)
# no slide down: -Fpar <= mu N ; no slide up: Fpar <= mu N
low=sp.solve(sp.Eq(-Fpar,mu*Nn),w**2)  # w^2 lower bound
up=sp.solve(sp.Eq(Fpar,mu*Nn),w**2)
print('low w^2=',sp.simplify(low[0]) if low else low)
print('up  w^2=',sp.simplify(up[0]) if up else up)
wl2=g*(sp.cos(th)-mu*sp.sin(th))/(s*sp.sin(th)*(sp.sin(th)+mu*sp.cos(th)))
wu2=g*(sp.cos(th)+mu*sp.sin(th))/(s*sp.sin(th)*(sp.sin(th)-mu*sp.cos(th)))
assert sp.simplify(low[0]-wl2)==0 and sp.simplify(up[0]-wu2)==0
# brute force check of classification
rng=np.random.default_rng(0)
for _ in range(20000):
    thv=rng.uniform(0.05,1.52); muv=rng.uniform(0.02,0.98); sv=1.0; gv=9.8; wv=rng.uniform(0,30)
    rv=sv*np.sin(thv)
    fp=rv*wv**2*np.sin(thv)-gv*np.cos(thv); nn=rv*wv**2*np.cos(thv)+gv*np.sin(thv)
    stay= abs(fp)<=muv*nn
    tl=np.tan(thv)
    lowok = True if tl>=1/muv else wv**2>=gv*(np.cos(thv)-muv*np.sin(thv))/(sv*np.sin(thv)*(np.sin(thv)+muv*np.cos(thv)))
    upok  = True if tl<=muv else wv**2<=gv*(np.cos(thv)+muv*np.sin(thv))/(sv*np.sin(thv)*(np.sin(thv)-muv*np.cos(thv)))
    if abs(abs(fp)-muv*nn)>1e-9: assert stay==(lowok and upok)
# (3) theta=45deg, fixed w: allowed s in [s1,s2]; s2/s1=9 -> mu
s1,s2=sp.symbols('s1 s2',positive=True)
th45={th:sp.pi/4}
S1=sp.solve(sp.Eq(w**2,wl2.subs(th45)),s)[0]; S2=sp.solve(sp.Eq(w**2,wu2.subs(th45)),s)[0]
ratio=sp.simplify(S2/S1); print('s2/s1=',ratio)
musol=sp.solve(sp.Eq(ratio,9),mu); print('mu=',musol)
wsol=sp.solve(sp.Eq(S1.subs(mu,sp.Rational(1,2)),s1),w)[0]; print('w=',sp.simplify(wsol))
# friction-free position s0
S0=sp.solve(sp.Eq(Fpar.subs(th45),0),s)[0]; print('s0/s1=',sp.simplify(S0/S1.subs(mu,sp.Rational(1,2))))
print('OK')
