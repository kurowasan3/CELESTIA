# N35-1: cylinder radius A (axis O at origin) with cylindrical cavity radius b, axis O'(d,0); uniform current density J (+z), total I.
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
A,b,d=1.0,0.35,0.45; I=1.0; J=I/(np.pi*(A**2-b**2))
# numeric superposition of filaments on fine polar-ish grid: H = sum dI/(2 pi r) phi-hat
n=1600; xs=np.linspace(-A,A,n); dx=xs[1]-xs[0]
X,Y=np.meshgrid(xs,xs); mask=(X**2+Y**2<=A**2)&((X-d)**2+Y**2>=b**2)
px,py=X[mask],Y[mask]; dI=J*dx*dx
def H(x,y):
    rx,ry=x-px,y-py; r2=rx**2+ry**2
    return np.array([np.sum(-ry/r2),np.sum(rx/r2)])*dI/(2*np.pi)
print('theory H inside cavity: (0, J d/2)=',(0,J*d/2))
for pt in [(d,0),(d+0.2,0.1),(d-0.25,-0.1),(d,0.3)]:
    print(pt,H(*pt))
# symbolic: H_full(r)=J r/2 (r<A) ; H_cav = (J/2) z x (r - r') = (J/2) z x d
Jj,dd,x,y=sp.symbols('J d x y')
Hf=sp.Matrix([-y,x])*Jj/2; Hh=-sp.Matrix([-y,x-dd])*Jj/2
print('symbolic cavity field:',sp.simplify(Hf+Hh).T)
# particle: q,m, B=mu0 H yhat uniform, start at O'(d,0,0) velocity v zhat. orbit in x-z plane; max x-deviation 2 m v/(qB)
q,m,Bc,v=1.0,1.0,2.0,1.0
def rhs(t,u):
    vel=u[3:]; F=q*np.cross(vel,[0,Bc,0])/m
    return np.concatenate([vel,F])
s=solve_ivp(rhs,[0,2*np.pi*m/(q*Bc)],[d,0,0,0,0,v],rtol=1e-10,atol=1e-12,dense_output=True)
print('min x =',s.y[0].min(),' expected d-2mv/qB=',d-2*m*v/(q*Bc),' max z=',s.y[2].max(),' expected mv/qB=',m*v/(q*Bc))
print('back to start after T:',s.y[:3,-1])
