# T20201-1: puck on frictionless table tied through hole to hanging mass M
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,M,g,A,r,a=sp.symbols('m M g A r a',positive=True)
U=2*m*A**2/r**2+M*g*r     # effective potential (KE of rotation + PE of hanging mass, const dropped)
r0=sp.solve(sp.diff(U,r),r)[0]; print('r0=',r0, ' check (4mA^2/Mg)^(1/3):',sp.simplify(r0**3-4*m*A**2/(M*g)))
print('U\'\'(r0)=',sp.simplify(sp.diff(U,r,2).subs(r,r0)*r0/(M*g)),'* Mg/r0')
Emin=sp.simplify(U.subs(r,r0)); print('Emin=',Emin, sp.simplify(Emin/(M*g*r0)))
wr2=sp.simplify(sp.diff(U,r,2).subs(r,r0)/(m+M)); wphi2=sp.simplify((2*A/r0**2)**2)
print('wr^2/wphi^2=',sp.simplify(wr2/wphi2)); print('wphi^2*r0/g=',sp.simplify(wphi2*r0/g))
# launch at r=a with tangential speed v=2 sqrt(ga), M=m
v=2*sp.sqrt(g*a); Av=a*v/2
E0=sp.Rational(1,2)*m*v**2+m*g*a
rr=sp.symbols('rr',positive=True)
print('turning pts:',sp.solve(sp.Eq(2*m*Av**2/rr**2+m*g*rr,E0),rr))
# tension formula T = M(g + rddot), rddot=(4mA^2/r^3 - Mg)/(m+M)
for R_ in [a,a*(1+sp.sqrt(3))]:
    rdd=(4*m*Av**2/R_**3-m*g)/(2*m); print('r=',R_,' rddot/g=',sp.nsimplify(sp.simplify(rdd/g)),float(rdd/g),' T/mg=',float(sp.simplify(m*(g+rdd)/(m*g))))
# speed at outer point
rmax=a*(1+sp.sqrt(3)); vmax_t=sp.simplify(2*Av/rmax); print('speed at rmax:',vmax_t, sp.simplify(vmax_t/sp.sqrt(g*a)))
# numeric simulation in polar coords with constraint: (m+M) r'' = m r phi'^2 - M g ; (r^2 phi')'=0
gv,av,mv,Mv=9.8,1.0,1.0,1.0
def rhs(t,y):
    R,Rd,ph,phd=y
    return [Rd,(mv*R*phd**2-Mv*gv)/(mv+Mv),phd,-2*Rd*phd/R]
y0=[av,0,0,2*np.sqrt(gv*av)/av]
s=solve_ivp(rhs,[0,10],y0,max_step=1e-3,rtol=1e-11,atol=1e-12)
print('numeric rmax=',s.y[0].max(),' pred',1+np.sqrt(3),' rmin=',s.y[0].min())
Rdd=(mv*s.y[0]*s.y[3]**2-Mv*gv)/(mv+Mv); T=Mv*(gv+Rdd); print('min tension/mg',T.min()/(mv*gv),' max',T.max()/(mv*gv))
# check: tension acting on puck: m(r''-r phi'^2) = -T
print('consistency',np.max(np.abs(mv*(Rdd-s.y[0]*s.y[3]**2)+T)))
# small oscillation period check, M=2m -> wr=wphi
Mv=2.0; Av_=0.7; r0v=(4*mv*Av_**2/(Mv*gv))**(1/3)
y0=[r0v*1.01,0,0,2*Av_/(r0v*1.01)**2]
s=solve_ivp(rhs,[0,30],y0,max_step=1e-3,rtol=1e-11,atol=1e-12,dense_output=True)
t=np.linspace(0,30,300001); R=s.sol(t)[0]
# radial period from successive maxima
from scipy.signal import argrelmax
idx=argrelmax(R)[0]; Tr=np.mean(np.diff(t[idx])); Tphi=2*np.pi/(2*Av_/r0v**2)
print('M=2m: Tr=',Tr,' Tphi=',Tphi)
