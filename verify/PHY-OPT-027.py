# J087-1: two convex lenses separated by d; magnification independence of object position (afocal), longitudinal mag, 4f relay
import sympy as sp, numpy as np
x,f1,f2,d=sp.symbols('x f1 f2 d',positive=True)
b1=f1*x/(x-f1); a2=d-b1; b2=sp.simplify(f2*a2/(a2-f2))
m=sp.simplify((-b1/x)*(-b2/a2))
print('m general =',sp.factor(m))
print('dm/dx numerator:',sp.factor(sp.numer(sp.together(sp.diff(m,x)))))
mA=sp.simplify(m.subs(d,f1+f2)); b2A=sp.simplify(b2.subs(d,f1+f2))
print('afocal m=',mA,' b2=',sp.expand(b2A))
print('db2/dx=',sp.simplify(sp.diff(b2A,x)))
f=sp.symbols('f',positive=True)
print('relay f1=f2: b2=',sp.simplify(b2A.subs({f1:f,f2:f})),' total distance x+2f+b2 =',sp.simplify(x+2*f+b2A.subs({f1:f,f2:f})))
# matrix-method cross check numeric (ray transfer)
def sysm(xv,f1v,f2v,dv):
    T=lambda L: np.array([[1,L],[0,1]]); Lm=lambda F: np.array([[1,0],[-1/F,1]])
    M=Lm(f2v)@T(dv)@Lm(f1v)@T(xv)
    # image distance z s.t. B element of T(z)M =0
    z=-M[0,1]/M[1,1]; Mt=T(z)@M; return z,Mt[0,0]
for xv in (0.05,0.12,0.3):
    z,mm=sysm(xv,0.1,0.25,0.35); print(xv,'z',z,'m',mm,' formula b2',float(b2A.subs({x:xv,f1:0.1,f2:0.25})),' m',float(mA.subs({f1:0.1,f2:0.25})))
