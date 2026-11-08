import numpy as np, sympy as sp
from scipy.integrate import solve_ivp, quad
g=9.8; vt=4.0
x=lambda t:(vt**2/g)*np.log(np.cosh(g*t/vt))
print('table',[round(x(t),2) for t in (2.0,3.0,4.0,5.0)])
# fit late-time line
ts=np.array([2.,3.,4.,5.]); xs=np.array([round(x(t),2) for t in ts])
p=np.polyfit(ts,xs,1); print('slope,intercept',p, ' vt^2/g',vt**2/g,' ratio',-p[1]/(vt**2/g),' ln2',np.log(2))
# general n lag integral
for n in (1,2):
    I=quad(lambda s:(1-s)/(1-s**n) if s<1 else 1/n,0,1)[0]; print('n',n,'lag coeff',I)
# ODE check for n=2 and n=1 lag
for n in (1,2):
    c_over_m=g/vt**n
    sol=solve_ivp(lambda t,y:[y[1],g-c_over_m*y[1]**n],[0,20],[0,0],rtol=1e-11,atol=1e-12)
    print('n',n,'lag ODE',vt*20-sol.y[0][-1],' formula',(vt**2/g)*(1 if n==1 else np.log(2)))
# relaxation time general n: tau = vt/(n g)
v,c,m,n=sp.symbols('v c m n',positive=True); vts=(m*g/c)**(1/n)
print('tau check', sp.simplify(1/(sp.diff(c*v**n/m,v).subs(v,vts)) - vts/(n*g)))
# kick u=0.20 for n=2 exact extra distance
u=0.20
sol=solve_ivp(lambda t,y:[y[1],g-(g/vt**2)*y[1]**2],[0,20],[0,vt+u],rtol=1e-12,atol=1e-12)
print('extra dist ODE',sol.y[0][-1]-vt*20,' exact',(vt**2/g)*np.log(1+u/(2*vt)),' linear',u*vt/(2*g))
print('2m: vt',vt*np.sqrt(2),' lag',2*(vt**2/g)*np.log(2))
