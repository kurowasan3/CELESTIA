# T20232-1: falling coil in radial field with counterweight; two runs eliminate BL and coil resistance
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,M,g,K,R,r,R1,R2,t=sp.symbols('m M g K R r R_1 R_2 t',positive=True)
# terminal: (M+m)g - Mg - K I = 0 ; K v = I (R+r)
I=m*g/K; v=I*(R+r)/K
print('I =',I,'  v =',sp.simplify(v))
v1=v.subs(R,R1); v2=v.subs(R,R2); V1=I*R1; V2=I*R2
m_expr=(V1/R1)**2*(R1-R2)/(g*(v1-v2)); print('m formula check:',sp.simplify(m_expr-m))
r_expr=(R1*v2-R2*v1)/(v1-v2); print('r formula check:',sp.simplify(r_expr-r))
print('power balance mgv=I^2(R+r):',sp.simplify(m*g*v-I**2*(R+r)))
# time constant
tau=(2*M+m)*(R+r)/K**2; print('tau*m g / v =',sp.simplify(tau*m*g/v))
# numbers
R1n,R2n,V1n,V2n,v1n,v2n,gn=100.,50.,0.200,0.100,3.00e-3,2.00e-3,9.80
In=V1n/R1n; Kn=In*(R1n-R2n)/(v1n-v2n); mn=Kn*In/gn; rn=(R1n*v2n-R2n*v1n)/(v1n-v2n)
print('I',In,' K',Kn,' m',mn,' r',rn)
# numerical ODE with those numbers (choose M=5 kg)
Mn=5.0
f=lambda tt,y:[(mn*gn-Kn**2*y[0]/(R1n+rn))/(2*Mn+mn)]
s=solve_ivp(f,[0,1],[0],rtol=1e-10,atol=1e-14,dense_output=True)
vt=mn*gn*(R1n+rn)/Kn**2; taun=(2*Mn+mn)*(R1n+rn)/Kn**2
print('terminal',s.y[0,-1],vt,' v(tau)/vt',s.sol(taun)[0]/vt,1-np.exp(-1))

# (5) new: v(t)=v1(1-exp(-t/tau)), tau=(2M+m)(R1+r)/K^2 ; M cancels in weight but not inertia
vf=sp.Function('v')
sol=sp.dsolve(sp.Eq((2*M+m)*vf(t).diff(t), m*g-K**2*vf(t)/(R1+r)), vf(t), ics={vf(0):0}).rhs
tau1s=(2*M+m)*(R1+r)/K**2
print('v(t) check:',sp.simplify(sol - v1*(1-sp.exp(-t/tau1s))))
print('tau = (2M+m) v1/(m g):',sp.simplify(tau1s-(2*M+m)*v1/(m*g)))
print('v1 independent of M:',sp.diff(v1,M)==0,'  tau depends on M: dtau/dM =',sp.diff(tau1s,M))
# inductance: L dI/dt + (R1+r) I = K v ; (2M+m) dv/dt = m g - K I
tauL_crit=(2*Mn+mn)*(R1n+rn)**2/Kn**2
print('L scale (2M+m)(R1+r)^2/K^2 =',tauL_crit,'H')
def runL(Lv):
    f=lambda tt,y:[(mn*gn-Kn*y[1])/(2*Mn+mn),(Kn*y[0]-(R1n+rn)*y[1])/Lv]
    tt=np.linspace(0,8*taun,4001)
    s=solve_ivp(f,[0,tt[-1]],[0,0],t_eval=tt,rtol=1e-10,atol=1e-14,method='LSODA')
    dev=np.max(np.abs(s.y[0]-vt*(1-np.exp(-tt/taun))))/vt
    return dev, s.y[0].max()/vt
for frac in [1e-3,1e-2,0.1,1.0,10.0]:
    dev,mx=runL(frac*tauL_crit)
    print(f'L/[(2M+m)(R1+r)^2/K^2]={frac:g}: max rel deviation from exponential={dev:.3e}, max v/v_terminal={mx:.4f}')
dev,_=runL(1e-3*tauL_crit); assert dev<5e-3
_,mx=runL(10*tauL_crit); assert mx>1.0   # overshoot when L large
