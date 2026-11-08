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
tau1=0.0153; print('2M+m = m g tau1/v1 =',mn*gn*tau1/v1n,' M=',(mn*gn*tau1/v1n-mn)/2)
# numerical ODE with those numbers (choose M=5 kg)
Mn=5.0
f=lambda tt,y:[(mn*gn-Kn**2*y[0]/(R1n+rn))/(2*Mn+mn)]
s=solve_ivp(f,[0,1],[0],rtol=1e-10,atol=1e-14,dense_output=True)
vt=mn*gn*(R1n+rn)/Kn**2; taun=(2*Mn+mn)*(R1n+rn)/Kn**2
print('terminal',s.y[0,-1],vt,' v(tau)/vt',s.sol(taun)[0]/vt,1-np.exp(-1))
