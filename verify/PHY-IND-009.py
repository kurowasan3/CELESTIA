# N44-1: inverse problem: coil with internal resistance; pickup coil isolates inductive emf
import numpy as np, sympy as sp
from scipy.integrate import quad
I0=2.0; Vflat=0.60; VQ=0.020; N1=2000; N2=80; l=0.50; mu0=4e-7*np.pi
k_up=I0/0.40; k_dn=-I0/0.10
r=Vflat/I0; M=VQ/k_up; L=M*N1/N2; S=L*l/(mu0*N1**2)
print('r=',r,' M=',M,' L=',L,' S=',S)
I=lambda t: k_up*t if t<0.4 else (I0 if t<0.8 else (I0+k_dn*(t-0.8) if t<0.9 else 0))
dI=lambda t: k_up if t<0.4 else (0 if t<0.8 else (k_dn if t<0.9 else 0))
V=lambda t: L*dI(t)+r*I(t)
for tt in [1e-9,0.4-1e-9,0.6,0.8+1e-9,0.9-1e-9]: print('V(%.3f)=%.4f'%(tt,V(tt)))
print('VQ during ramp-down magnitude',abs(M*k_dn))
Wdn=quad(lambda t:V(t)*I(t),0.8,0.9,points=[0.8,0.9])[0]
Wtot=sum(quad(lambda t:V(t)*I(t),a_,b_)[0] for a_,b_ in [(0,0.4),(0.4,0.8),(0.8,0.9)])
print('W(0.8-0.9)=',Wdn,'  W total=',Wtot,' Joule total=',sum(quad(lambda t:r*I(t)**2,a_,b_)[0] for a_,b_ in [(0,0.4),(0.4,0.8),(0.8,0.9)]))
# symbolic L/M=N1/N2 from solenoid formulas
mu,n1,n2,SS,ll=sp.symbols('mu0 N1 N2 S l',positive=True)
print('L/M =',sp.simplify((mu*n1**2*SS/ll)/(mu*n1*n2*SS/ll)))
