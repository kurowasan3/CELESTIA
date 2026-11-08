# K20261-1: neutral ball from free wedge merges with gyrating charge
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,M,g,h,v,u,V,R=sp.symbols('m M g h v u V R',positive=True)
s=sp.solve([sp.Eq(m*v-M*u,0),sp.Eq(m*v**2/2+M*u**2/2,m*g*h)],[v,u],dict=True)[0]
print('(1) v=',sp.simplify(s[v]))
# merged velocity (x, right +): A +v, C -V
vF=(m*v-m*V)/(2*m); print('(2) vF=',sp.simplify(vF))
q,B=sp.symbols('q B',positive=True)
RF=sp.simplify(2*m*sp.Abs(vF)/(q*B)).subs(q*B,m*V/R); print('    R_F=',sp.simplify(2*m*sp.Abs(vF)/(q*B)),'with qB=mV/R ->',sp.simplify(RF))
for vv in (V/2,3*V/2):
    hh=sp.solve(sp.Eq(s[v],vv),h)[0]; print('(3) v=',vv,' h=',sp.simplify(hh))
# (4) energy fraction for v=V/2
KEb=m*(V/2)**2/2+m*V**2/2; KEa=2*m*(V/4)**2/2; print('(4) KE after/before =',sp.simplify(KEa/KEb))
# numeric: Lorentz force integration of F, units m=1,q=1,B=1,V=1 -> R=1, O at (0,1), P at origin
def run(vA):
    vx=(vA-1.0)/2; mF=2.0; qF=1.0; Bz=1.0
    f=lambda t,y:[y[2],y[3],qF*(y[3]*Bz)/mF,qF*(-y[2]*Bz)/mF]
    sol=solve_ivp(f,[0,4*np.pi],[0,0,vx,0],rtol=1e-11,atol=1e-12,dense_output=True)
    ts=np.linspace(0,4*np.pi,200001); Y=sol.sol(ts)
    dO=np.hypot(Y[0],Y[1]-1.0); i=np.argmin(dO)
    return Y, ts[i], dO[i]
Y,tO,dmin=run(0.5)
print('    v=V/2: closest approach to O =',dmin,' at t=',tO,' (2πR/V=',2*np.pi,')', ' ymin,ymax=',Y[1].min(),Y[1].max())
Y2,_,d2=run(1.5); print('    v=3V/2: y range',Y2[1].min(),Y2[1].max(),' x range',Y2[0].min(),Y2[0].max())
# check C's original sense: q>0, B up, velocity -x at P(0,0) -> force direction
vC=np.array([-1,0,0]); Bv=np.array([0,0,1]); print('    force on C at P:',np.cross(vC,Bv),'(+y: toward O, ok)')
# rotation sense of F (v=V/2): angular momentum about its centre (0,0.5)
x,y=Y[0][1000],Y[1][1000]; 
