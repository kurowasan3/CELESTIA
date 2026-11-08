# T20093-1: piston on coexisting liquid+vapor. Check (a) energy balance -> C dT = q dn_cond
# (b) small-oscillation period under adiabatic condition by direct ODE integration.
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
# symbolic first-order check
C,q,p,dv,A,y,al,M=sp.symbols('C q p Delta_v A y alpha M',positive=True)
dT=sp.symbols('dT')
dn_v=-A*y/dv           # vapor mol change when piston pushed down by y
uL=q-p*dv               # internal energy of vaporization per mol
sol=sp.solve(sp.Eq(C*dT+dn_v*uL, p*A*y),dT)[0]
print('dT =',sp.simplify(sol))
k=sp.simplify(A*al*sol/y); print('k_eff =',k, ' period=',sp.simplify(2*sp.pi*sp.sqrt(M/k)))
# isothermal heat out
print('Q_out(isothermal)=',sp.simplify(-(dn_v*uL - p*A*y)))
# numeric ODE
Cn,qn,p0,al_n,A_n,vg,vl,M_n,g=500.0,4.4e4,4.2e3,250.0,1e-2,0.6,1.8e-5,None,9.8
M_n=p0*A_n/g; dvn=vg-vl; uLn=qn-p0*dvn; T0=303.0
def rhs(t,s):
    z,zd,T=s  # z: piston height displacement (up +)
    pr=p0+al_n*(T-T0)
    zdd=(pr*A_n-M_n*g)/M_n
    dV=A_n*zd; dn=dV/dvn
    dTdt=-(uLn+pr*dvn)*dn/Cn
    return [zd,zdd,dTdt]
z0=-1e-4
T_init=T0+qn*A_n*(-z0)/(Cn*dvn)
s=solve_ivp(rhs,[0,200],[z0,0,T_init],max_step=0.01,rtol=1e-10,atol=1e-14)
z=s.y[0]; t=s.t
cross=[t[i] for i in range(1,len(t)) if z[i-1]<0<=z[i]]
Tnum=np.mean(np.diff(cross))
Tform=2*np.pi*np.sqrt(M_n*Cn*dvn/(al_n*qn*A_n**2))
print('period numeric',Tnum,'formula',Tform)
