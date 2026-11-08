# N20-1: charged isolated plate M between L (x=-a) and N (x=+a); L grounded, N at potential V0 (V0=0 for parts 1-4)
import sympy as sp, numpy as np
a,x,q,V0,C0,k,m,v=sp.symbols('a x q V0 C0 k m v',real=True)
eS=C0*a                          # eps*S
C1=eS/(a+x); C2=eS/(a-x)
VM=sp.symbols('VM')
VMs=sp.solve(sp.Eq(C1*(VM-0)+C2*(VM-V0),q),VM)[0]
QL=sp.simplify(C1*(0-VMs)); QN=sp.simplify(C2*(V0-VMs))
print('VM',sp.factor(VMs)); print('QL',sp.factor(QL)); print('QN',sp.factor(QN))
# current into N when M moves at v (dx/dt=v)
print('dQN/dx', sp.simplify(sp.diff(QN,x)))
# force on M (+x): faces: right face field E_R=(VM-V0)/(a-x), left face E_L=VM/(a+x); pressure eps E^2/2 outward
ER=(VMs-V0)/(a-x); EL=VMs/(a+x)
F=sp.simplify(eS/2*(ER**2-EL**2)/1)  # eps*S/2 * E^2 , eps S = eS (S cancels: eps*E^2/2 * S)
print('F',sp.expand(sp.simplify(F)))
# energy-method check (V0=0): U = q^2/(2(C1+C2)), F=-dU/dx
U=q**2/(2*(C1+C2)); print('F energy (V0=0)', sp.simplify(-sp.diff(U,x)), ' vs ', sp.simplify(F.subs(V0,0)))
# energy-method with battery: Fx = -d/dx[ q^2/(2(C1+C2)) ... ] -- use generalized: F = -dU/dx|_q + V0 dQN/dx ... check
Ue=sp.Rational(1,2)*(C1*VMs**2+C2*(V0-VMs)**2)
Fe=sp.simplify(-sp.diff(Ue,x)+V0*sp.diff(QN,x))
print('F energy w/ battery minus F', sp.simplify(Fe-F))
# numeric oscillation check: m=1,k=3,q=1,a=1,C0=1 -> keff = k - q^2/(2a^2C0)=2.5
from scipy.integrate import solve_ivp
Ff=sp.lambdify((x,),F.subs({q:1,a:1,C0:1,V0:0}))
f=lambda t,y:[y[1],-3*y[0]+Ff(y[0])]
s=solve_ivp(f,[0,30],[1e-3,0],max_step=1e-3,rtol=1e-11,atol=1e-14)
yy=s.y[0]; idx=np.where((yy[:-1]<0)&(yy[1:]>=0))[0]; T=np.mean(np.diff(s.t[idx]))
print('omega num',2*np.pi/T,'theory',np.sqrt(2.5))
