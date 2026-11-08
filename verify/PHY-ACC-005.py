import sympy as sp, numpy as np
from scipy.optimize import minimize_scalar, brentq
V0,R,wb,ws,Ib=10.0,20.0,1.0e4,2.0e4,0.50
L2=V0/(wb*Ib); C=1/(wb**2*L2); L1=L2*wb**2/(ws**2-wb**2)
print('L2,C,L1',L2,C,L1)
def Z(w): return R+1j*w*L1+1/(1/(1j*w*L2)+1j*w*C)
w=np.linspace(100,1e5,2000001)
I=np.abs(V0/Z(w))
print('zero of I near', w[np.argmin(I)], ' peak near', w[np.argmax(I)], 'peak',I.max())
# circulating current at wb: tank voltage
w0=wb*(1+1e-9); Itot=V0/Z(w0); Vt=Itot/(1/(1j*w0*L2)+1j*w0*C)
print('|I_L2| at wb', abs(Vt/(1j*w0*L2)))
It=V0/Z(ws); Vt=It/(1/(1j*ws*L2)+1j*ws*C)
print('at ws: I',abs(It),' V_L1',abs(It*1j*ws*L1),' Vtank',abs(Vt),' IC',abs(Vt*1j*ws*C),' IL2',abs(Vt/(1j*ws*L2)))
X=np.imag(Z(w)); lead=w[X<0]
print('current leads for w in',lead.min(),lead.max())
print('low/high limits', abs(V0/Z(1.0)), abs(V0/Z(1e9)))
# symbolic series resonance
l1,l2,c,om=sp.symbols('L1 L2 C omega',positive=True)
Xs=om*l1+om*l2/(1-om**2*l2*c)
print(sp.solve(sp.Eq(sp.together(Xs).as_numer_denom()[0],0),om))
