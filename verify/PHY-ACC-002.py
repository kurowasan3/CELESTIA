# N48-1: delta network of R, L, C boxes; identify from step responses
import sympy as sp, itertools, numpy as np
from scipy.integrate import solve_ivp
s=sp.symbols('s',positive=True)
E,r=12,4
Rv,Lv,Cv=8,sp.Rational(8,100),sp.Rational(5,100000)
Z={'R':Rv,'L':s*Lv,'C':1/(s*Cv)}
def par(a,b): return a*b/(a+b)
res={}
for perm in itertools.permutations('RLC'):
    X,Y,Zz=[Z[e] for e in perm]   # X:AB, Y:BC, Z:CA
    pairs={'AB':par(X,Y+Zz),'BC':par(Y,X+Zz),'CA':par(Zz,X+Y)}
    tab=[]
    for k_,Zp in pairs.items():
        I=E/(s*(r+Zp))
        i0=sp.limit(s*I,s,sp.oo); iinf=sp.limit(s*I,s,0)
        tab.append((k_,i0,iinf))
    res[perm]=tab
    print(perm,[(a,float(b),float(c)) for a,b,c in tab])
target=[('AB',1,1),('BC',1,3),('CA',3,1)]
match=[p for p,t in res.items() if [(a,float(b),float(c)) for a,b,c in t]==[(a,float(b),float(c)) for a,b,c in target]]
print('matching assignments:',match)
assert match==[('R','L','C')]
# R from E/(r+R)=1
print('R =',sp.solve(sp.Eq(sp.Integer(E)/(r+sp.Symbol('R')),1))[0])
# AC between A,B with r: |r+Z_AB(jw)| min at w0
w=sp.symbols('w',positive=True); j=sp.I
ZAB=par(Rv,j*w*Lv+1/(j*w*Cv))
mag=sp.lambdify(w,sp.Abs(r+ZAB))
ws=np.linspace(100,1000,90001); m=np.array([mag(x) for x in ws]); print('min |r+Z| at w=',ws[m.argmin()],' value',m.min())
# opening switch in BC steady state: coil current 3 A into loop L-R-C (C uncharged)
L_,C_,R_=float(Lv),float(Cv),float(Rv)
def f(t,y):
    i,vc,H=y; return [(-R_*i-vc)/L_, i/C_, R_*i*i]
sol=solve_ivp(f,[0,0.5],[3,0,0],rtol=1e-10,atol=1e-12,max_step=1e-5)
print('V_L just after =',R_*3,' heat =',sol.y[2,-1],' 1/2 L I^2 =',0.5*L_*9)
print('L from heat 0.36 J:',2*0.36/9,' C from w0=500:',1/(500**2*0.08))
