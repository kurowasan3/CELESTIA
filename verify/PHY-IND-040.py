# Q171-1 three lamps, coils L and 2L, after opening the switch
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
L,r,E=sp.symbols('L r E',positive=True); I0=E/r
Q1,Q2=sp.symbols('Q1 Q2')
# L1 di1/dt=-r(2i1-i2), L2 di2/dt=r(i1-2i2); L1=L,L2=2L; i1(0)=2I0,i2(0)=I0
s=sp.solve([sp.Eq(-L*2*I0,-r*(2*Q1-Q2)),sp.Eq(-2*L*I0,r*(Q1-2*Q2))],[Q1,Q2]); print(s, 'lamp2:',sp.simplify(s[Q1]-s[Q2]))
print('total heat', sp.Rational(1,2)*L*(2*I0)**2+sp.Rational(1,2)*2*L*I0**2)
# numeric
Lv,rv,Ev=1.0,1.0,1.0; I=Ev/rv
f=lambda t,y:[-rv*(2*y[0]-y[1])/Lv, rv*(y[0]-2*y[1])/(2*Lv)]
t=np.linspace(0,40,400001)
so=solve_ivp(f,[0,40],[2*I,I],t_eval=t,rtol=1e-11,atol=1e-13)
i1,i2=so.y; i_l2=i1-i2
from numpy import trapezoid as tr
print('Q lamp1',tr(i1,t),'Q lamp3',tr(i2,t),'Q lamp2',tr(i_l2,t))
print('heat total',tr(rv*(i1**2+i_l2**2+i2**2),t),'theory',3*Lv*I**2)
k=np.where(np.sign(i_l2[:-1])!=np.sign(i_l2[1:]))[0]; print('lamp2 sign change at t=',t[k][:3])
print('lamp heats', tr(rv*i1**2,t),tr(rv*i_l2**2,t),tr(rv*i2**2,t))
# eigen check
M=sp.Matrix([[-2,1],[sp.Rational(1,2),-1]]); print(M.eigenvals())
# reversal time analytic
lam=sorted([complex(x).real for x in M.eigenvals()]); print(lam)
