# N67-1: d + 3He -> 4He + p, proton at 90 deg determines Q and B(3He); mirror nuclei Coulomb
import sympy as sp, numpy as np
Kd,Kp,Ka,Q=sp.symbols('K_d K_p K_a Q')
# masses ~ mass numbers: d=2, p=1, a=4 ; momentum^2 = 2 m K
eqs=[sp.Eq(2*2*Kd+2*1*Kp,2*4*Ka), sp.Eq(Kd+Q,Kp+Ka)]
s=sp.solve(eqs,[Ka,Q],dict=True)[0]; print('Q =',sp.simplify(s[Q]))
Qv=float(s[Q].subs({Kd:2.00,Kp:15.49})); print('Q=',Qv)
Bd,Ba,Bt=2.22,28.30,8.48
BHe3=Ba-Bd-Qv; print('B(3He)=',BHe3,' diff',Bt-BHe3)
# vector check: build momenta
pd=np.array([np.sqrt(2*2*2.00),0]); pp=np.array([0,np.sqrt(2*1*15.49)]); pa=pd-pp
Kav=pa@pa/(2*4); print('K_alpha',Kav,' energy check Q', Kav+15.49-2.00)
print('alpha angle below axis (deg)',np.degrees(np.arctan2(-pa[1],pa[0])))
k0,e=9.0e9,1.6e-19
dE=(Bt-BHe3)*1.6e-13; print('r=',k0*e**2/dE)
