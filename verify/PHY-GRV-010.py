# K20231-1 verification: point mass + uniform dust sphere, radial oscillation & apsidal angle
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
G,M,R,rho,r,m=sp.symbols('G M R rho r m',positive=True)
k=sp.Rational(4,3)*sp.pi*G*rho
f=G*M/r**2+k*r            # attraction per unit mass
Om2=sp.simplify(f.subs(r,R)/R)
h=R**2*sp.sqrt(Om2)
# radial eq: r''=h^2/r^3 - f(r)
g=h**2/r**3-f
wr2=sp.simplify(-sp.diff(g,r).subs(r,R))
print('Omega^2=',sp.simplify(Om2),' wr^2=',sp.expand(wr2))
assert sp.simplify(wr2-(G*M/R**3+4*k))==0
# inverse: Phi=2pi-delta
d=sp.symbols('delta',positive=True); x=sp.symbols('x',positive=True)
q=1-d/(2*sp.pi)
xs=sp.solve(sp.Eq((1+x)/(1+4*x),q**2),x)[0]
print('x=',sp.simplify(xs),' series:',sp.series(xs,d,0,2))
rho_sol=sp.simplify(3*M*xs/(4*sp.pi*R**3))
print('rho=',rho_sol,' approx',sp.series(rho_sol,d,0,2))
# numeric check of apsidal angle
Gn,Mn,Rn=1.0,1.0,1.0; rhon=0.05
kn=4/3*np.pi*Gn*rhon
Omn=np.sqrt(Gn*Mn/Rn**3+kn); wrn=np.sqrt(Gn*Mn/Rn**3+4*kn)
u0=1e-3
def rhs(t,y):
    x,yy,vx,vy=y; rr=np.hypot(x,yy); a=-(Gn*Mn/rr**3+kn)
    return [vx,vy,a*x,a*yy]
sol=solve_ivp(rhs,[0,60],[Rn,0,u0,Omn*Rn],rtol=1e-11,atol=1e-13,dense_output=True,max_step=0.01)
t=sol.t; rr=np.hypot(sol.y[0],sol.y[1]); ph=np.unwrap(np.arctan2(sol.y[1],sol.y[0]))
idx=[i for i in range(1,len(t)-1) if rr[i]<rr[i-1] and rr[i]<rr[i+1]]
dphi=np.diff(ph[idx])
print('numeric apsidal angle',dphi[:4],' theory',2*np.pi*Omn/wrn)
print('amplitude numeric',(rr.max()-rr.min())/2,' theory',u0/wrn)
assert abs(dphi.mean()-2*np.pi*Omn/wrn)<2e-3
# inverse check: delta from numeric -> rho
delta=2*np.pi-dphi.mean()
xx=float(xs.subs(d,delta)); print('rho recovered',3*Mn*xx/(4*np.pi*Rn**3),' true',rhon)
