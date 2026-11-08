# K20022-1 verification: rod on V-shaped rails in 1/x field, rod resistance per length lambda
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
mu,I,k,lam,m,v0,x0,x,x1=sp.symbols('mu0 I k lambda m v0 x0 x x1',positive=True)
B=mu*I/(2*sp.pi*x)
Phi=sp.integrate(B*k*sp.Symbol('s'),(sp.Symbol('s'),x1,x)) if False else sp.integrate((mu*I/(2*sp.pi*sp.Symbol('s')))*k*sp.Symbol('s'),(sp.Symbol('s'),x1,x))
print('flux:',sp.simplify(Phi),' dPhi/dx:',sp.simplify(sp.diff(Phi,x)))
K0=mu*I*k/(2*sp.pi); G=K0**2/(lam*k)
xf=x0*sp.exp(m*v0/G); print('x_f =',sp.simplify(xf))
print('Q =',sp.simplify(K0/(lam*k)*sp.log(xf/x0)),' m v0/K0 =',sp.simplify(m*v0/K0))
print('x at v0/2:',sp.simplify(x0*sp.exp(m*v0/(2*G))))
# numeric ODE
vals={mu:4e-7*np.pi,I:2*np.pi/(4e-7*np.pi)/0.5*0.3,k:0.5,lam:0.4,m:0.05,v0:2.0,x0:0.1}  # dimensionless-ish: K0=0.3
K0n=float(K0.subs(vals)); Gn=float(G.subs(vals)); mn=0.05
def rhs(t,y):
    X,V,Q,H=y
    i=K0n*V/(vals[lam]*vals[k]*X)
    return [V,-i*K0n/mn,i,i*i*vals[lam]*vals[k]*X]
sol=solve_ivp(rhs,[0,2000],[0.1,2.0,0,0],rtol=1e-10,atol=1e-14)
Xe,Ve,Qe,He=sol.y[:,-1]
print('Gamma=',Gn,' num x_end',Xe,' formula xf',float(xf.subs(vals)),' v_end',Ve)
print('num Q',Qe,' formula',float((m*v0/K0).subs(vals)),' heat',He,' 1/2mv0^2',0.5*mn*4)
# check v(x) relation at an intermediate time
j=len(sol.t)//3; Xj,Vj=sol.y[0,j],sol.y[1,j]
print('v(x) check',Vj, 2.0-Gn/mn*np.log(Xj/0.1))
