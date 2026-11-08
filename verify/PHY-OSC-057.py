# PHY-OSC-057 verification: light right-angle frame with ball at vertex resting on two smooth pegs
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp

d,g,m,M,phi,phi0=sp.symbols('d g m M phi phi0',positive=True)
# geometry: O origin, P1(-d/2,0), P2(d/2,0), C on circle radius d/2 (Thales)
C=sp.Matrix([d/2*sp.sin(phi),-d/2*sp.cos(phi)])
P1=sp.Matrix([-d/2,0]); P2=sp.Matrix([d/2,0])
u1=(P1-C); u2=(P2-C)
print('(1) CP1.CP2 =',sp.simplify(u1.dot(u2)))  # 0 -> right angle at C
u1n=u1/sp.sqrt(sp.simplify(u1.dot(u1))); u2n=u2/sp.sqrt(sp.simplify(u2.dot(u2)))
# (2) fixed pegs: energy -> v^2 ; G radial
v2=g*d*(sp.cos(phi)-sp.cos(phi0))
G=m*v2/(d/2)+m*g*sp.cos(phi)
print('(2) G =',sp.simplify(G))
# decomposition: G*(O-C)/|OC| = N1*u2n + N2*u1n
eOC=-C/(d/2)
N1=sp.simplify(G*eOC.dot(u2n)); N2=sp.simplify(G*eOC.dot(u1n))
for val in [0.3,-0.5,1.0]:
    s={phi:val,d:1}
    print('phi',val,' N1/G',float((eOC.dot(u2n)).subs(s)),' cos(pi/4+phi/2)=',np.cos(np.pi/4+val/2),
          ' N2/G',float((eOC.dot(u1n)).subs(s)),' cos(pi/4-phi/2)=',np.cos(np.pi/4-val/2))
# torque check about C: N1 acts at P1 along u2n, N2 acts at P2 along u1n; resultant through C?
def cross(a,b): return a[0]*b[1]-a[1]*b[0]
tq=sp.simplify(cross(P1-C,u2n)*eOC.dot(u2n)+cross(P2-C,u1n)*eOC.dot(u1n))
print('torque of peg forces about C (should be 0):',sp.simplify(tq))

# (4) cart of mass M: Lagrangian with X, phi
t=sp.symbols('t'); X=sp.Function('X')(t); f=sp.Function('f')(t)
xb=X+d/2*sp.sin(f); yb=-d/2*sp.cos(f)
T=sp.Rational(1,2)*M*X.diff(t)**2+sp.Rational(1,2)*m*(xb.diff(t)**2+yb.diff(t)**2)
V=m*g*yb; L=T-V
eqs=[sp.diff(L.diff(q.diff(t)),t)-L.diff(q) for q in (X,f)]
a1,a2=sp.symbols('a1 a2')
sol=sp.solve([e.subs({X.diff(t,2):a1,f.diff(t,2):a2}) for e in eqs],[a1,a2],dict=True)[0]
# small oscillation omega^2
fd=sp.symbols('fd')
acc=sp.simplify(sol[a2].subs({f.diff(t):0}).subs(f,phi))
w2=sp.simplify(-sp.diff(acc,phi).subs(phi,0))
print('(4) omega^2 =',w2,' expected 2(M+m)g/(M d):',sp.simplify(w2-2*(M+m)*g/(M*d)))
# numeric: d=1,g=9.8,m=1,M=2, phi0=1.2 (large) -> check max cart speed and peg forces positive
num={d:1.0,g:9.8,m:1.0,M:2.0}
fa1=sp.lambdify((f,fd),sol[a1].subs(f.diff(t),fd).subs(num))
fa2=sp.lambdify((f,fd),sol[a2].subs(f.diff(t),fd).subs(num))
def rhs(tt,y): return [y[1],fa1(y[2],y[3]),y[3],fa2(y[2],y[3])]
P0=1.2
s=solve_ivp(rhs,[0,5],[0,0,P0,0],max_step=1e-3,rtol=1e-10,atol=1e-12)
Vmax=np.max(np.abs(s.y[1]))
Vth=np.sqrt(1**2*9.8*1*(1-np.cos(P0))/(2*3))
print('(4) Vmax num',Vmax,' formula m*sqrt(g d(1-cos)/(M(M+m)))',Vth)
# check required frame force on ball is radial and N1,N2>0 along trajectory
dd=1.0; ok=True; maxdev=0
for k in range(0,len(s.t),200):
    Xv,Vv,fv,wv=s.y[:,k]
    ax=fa1(fv,wv); al=fa2(fv,wv)
    # ball accel
    bx=ax+dd/2*(np.cos(fv)*al-np.sin(fv)*wv**2)
    by=dd/2*(np.sin(fv)*al+np.cos(fv)*wv**2)
    F=np.array([1.0*bx,1.0*by+9.8])  # frame force on ball = m a - m g
    Cv=np.array([dd/2*np.sin(fv),-dd/2*np.cos(fv)])
    eoc=-Cv/np.linalg.norm(Cv)
    maxdev=max(maxdev,abs(F[0]*eoc[1]-F[1]*eoc[0]))
    p1=np.array([-dd/2,0])-Cv; p2=np.array([dd/2,0])-Cv
    U=np.column_stack([p2/np.linalg.norm(p2),p1/np.linalg.norm(p1)])
    n1,n2=np.linalg.solve(U,F)
    if n1<=0 or n2<=0: ok=False
print('frame force radial (max cross):',maxdev,' N1,N2>0 throughout:',ok)
# small amplitude period numeric
s2=solve_ivp(rhs,[0,6],[0,0,0.01,0],max_step=1e-3,rtol=1e-10,atol=1e-12,dense_output=True)
ff=s2.y[2]; idx=np.where((ff[:-1]>0)&(ff[1:]<=0))[0]
Tn=s2.t[idx[2]]-s2.t[idx[1]]
print('period num',Tn,' formula 2pi sqrt(M d/(2(M+m)g))',2*np.pi*np.sqrt(2*1/(2*3*9.8)))
