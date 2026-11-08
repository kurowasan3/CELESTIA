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

# ---- (4) rewritten: AC source V0 between terminals 1,2 ----
# (a) |I0| = V0/|r+Z12| <= V0/r for all w, equality only at w0=1/sqrt(LC), independent of R
Lsym,Csym,Rsym,rsym,V0=sp.symbols('L C R r V0',positive=True)
X=w*Lsym-1/(w*Csym)
Z12s=Rsym*(j*X)/(Rsym+j*X)
re12=sp.simplify(sp.re(sp.expand_complex(Z12s)))
print('Re Z12 =',sp.factor(re12),' (>=0 -> |r+Z12|>=r)')
assert sp.simplify(re12-Rsym*X**2/(Rsym**2+X**2))==0
w0=1/sp.sqrt(Lsym*Csym)
print('Z12 at w0:',sp.simplify(Z12s.subs(w,w0)))
assert sp.simplify(Z12s.subs(w,w0))==0
rng=np.random.default_rng(1)
for _ in range(2000):
    Rv_,rv_,Lv_,Cv_=rng.uniform(0.5,50),rng.uniform(0.5,50),rng.uniform(1e-3,1),rng.uniform(1e-6,1e-3)
    ww=rng.uniform(1,1e5)
    X_=ww*Lv_-1/(ww*Cv_); Z_=Rv_*1j*X_/(Rv_+1j*X_)
    assert abs(1/(rv_+Z_))<=1/rv_+1e-12
# (b) V_Z amplitude at resonance and energy constancy, time-domain check with the problem's R and r
Lf,Cf,Rf,rf,V0f=0.08,5e-5,8.0,4.0,1.0
w0f=1/np.sqrt(Lf*Cf)
def rhs(t,y):
    iL,vC=y            # LC-branch current, capacitor voltage ; node voltage v12 = R*(i - iL)
    # source loop: V0 sin = r*i + v12 ; branch: v12 = L diL/dt + vC
    i=(V0f*np.sin(w0f*t)+Rf*iL/1.0)/(rf+Rf)   # from V0 sin = r i + R(i-iL)
    v12=Rf*(i-iL)
    return [(v12-vC)/Lf, iL/Cf]
T0=2*np.pi/w0f
sol=solve_ivp(rhs,[0,400*T0],[0,0],max_step=T0/200,rtol=1e-9,atol=1e-12,dense_output=True)
tt=np.linspace(390*T0,400*T0,20001); y=sol.sol(tt)
iL,vC=y; i_src=(V0f*np.sin(w0f*tt)+Rf*iL)/(rf+Rf)
print('steady: I0 num',np.abs(i_src).max(),' V0/r',V0f/rf)
print('steady: V_Z num',np.abs(vC).max(),' formula',V0f/rf*np.sqrt(Lf/Cf))
U=0.5*Lf*iL**2+0.5*Cf*vC**2
print('stored energy min/max',U.min(),U.max(),' 1/2 L (V0/r)^2',0.5*Lf*(V0f/rf)**2)
assert abs(np.abs(vC).max()-V0f/rf*np.sqrt(Lf/Cf))<0.02*V0f/rf*np.sqrt(Lf/Cf)
assert abs(np.abs(i_src).max()-V0f/rf)<0.01*V0f/rf
assert (U.max()-U.min())<0.03*U.mean()
print('condition V_Z>V0 <=> L>r^2 C:', sp.solve(sp.Eq(V0/rsym*sp.sqrt(Lsym/Csym),V0),Lsym))
print('ACC-002 (4) OK')
