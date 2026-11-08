# S0201-1: C(T)=aT+bT^3 low-temperature metal sample
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
a,b,T,T1,T2=sp.symbols('a b T T1 T2',positive=True)
C=a*T+b*T**3
U=sp.integrate(C,(T,0,T))
print('U(T)=',U)
Q=sp.simplify(U.subs(T,T2)-U.subs(T,T1))
Qmean=sp.simplify(C.subs(T,(T1+T2)/2)*(T2-T1))
print('Q exact - Q(mean-T):',sp.factor(sp.simplify(Q-Qmean)))  # expect b(T2-T1)^3 (T1+T2)/4? check
# data generation
an,bn=1.0e-2,6.0e-4
Qp=1.00e3*(1.00e-3)**2*1.00
Tm=[1.50,2.50,3.50,4.50]
dT=[Qp/(an*t+bn*t**3) for t in Tm]
print('pulse Q',Qp,' dT(mK)',[round(x*1e3,1) for x in dT])
dTr=[round(x*1e3,1)/1e3 for x in dT]
Cm=[Qp/x for x in dTr]
X=np.array(Tm)**2; Y=np.array(Cm)/np.array(Tm)
p=np.polyfit(X,Y,1); print('fit from rounded data: b=',p[0],' a=',p[1])
print('C/T table:',[round(y*1e3,2) for y in Y],' T^2:',list(X))
# (2) 2.00->6.00 K
Qex=float(Q.subs({a:an,b:bn,T1:2,T2:6})); Qm=float(Qmean.subs({a:an,b:bn,T1:2,T2:6}))
print('Q exact',Qex,' Q meanT',Qm,' rel err',(Qm-Qex)/Qex)
print('electron part exact', an/2*(36-4),' meanT',an*4*4)
print('lattice part exact', bn/4*(6**4-2**4),' meanT',bn*64*4)
# (3) heater 1.00 mW 300 s from 2.00 K
P=1.00e-3
sol=solve_ivp(lambda t,y:[P/(an*y[0]+bn*y[0]**3)],[0,300],[2.0],rtol=1e-12,atol=1e-12)
x=sp.symbols('x',positive=True)
Uf=float(U.subs({a:an,b:bn,T:2}))+P*300
xs=sp.solve(sp.Eq(an/2*x+bn/4*x**2,Uf),x); print('T3 analytic',[sp.sqrt(v).evalf() for v in xs if v>0],' numeric',sol.y[0][-1])
print('naive with C(2.00K):', 2+P*300/(an*2+bn*8))
# (4) contact
Uc=float(U.subs({a:an,b:bn,T:2})+U.subs({a:an,b:bn,T:6}))/2
xs=sp.solve(sp.Eq(an/2*x+bn/4*x**2,Uc),x); Tf=[sp.sqrt(v).evalf() for v in xs if v>0][0]; print('Tf',Tf)
# numeric contact with conductance G
G=0.05
def rhs(t,y):
    q=G*(y[1]-y[0]); return [q/(an*y[0]+bn*y[0]**3), -q/(an*y[1]+bn*y[1]**3)]
sol=solve_ivp(rhs,[0,2000],[2.0,6.0],rtol=1e-11,atol=1e-12); print(' numeric contact final',sol.y[:,-1])
print('b=0 case Tf=',np.sqrt((4+36)/2))
# crossover temperature
print('T* = sqrt(a/b)=',np.sqrt(an/bn))
