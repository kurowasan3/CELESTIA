# Q108-1: hose roll unrolling up a slope (energy conservation with variable mass)
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp, quad
x,L,v0,g,b,lam=sp.symbols('x L v0 g beta lambda',positive=True)
u2=(L*v0**2-sp.Rational(4,3)*g*sp.sin(b)*(L*x-x**2/2))/(L-x)
E=sp.Rational(3,4)*lam*(L-x)*u2+lam*g*sp.sin(b)*(L*x-x**2/2)
print('energy const check:',sp.simplify(E-sp.Rational(3,4)*lam*L*v0**2))
xs=sp.solve(sp.Eq(L*v0**2-sp.Rational(4,3)*g*sp.sin(b)*(L*x-x**2/2),0),x); print('stop roots',xs)
# horizontal: time and force
u0=v0/sp.sqrt(1-x/L)
T=sp.integrate(1/u0,(x,0,L)); print('T horizontal',sp.simplify(T))
p=lam*(L-x)*u0               # momentum of whole hose (laid part at rest)
F=sp.simplify(sp.diff(p,x)*u0)  # dp/dt
print('F(x)',sp.factor(F))
# numeric ODE from energy: d/dt E=0 -> acceleration of roll; integrate
Lv,gv,bv,v0v=20.0,9.8,np.radians(10),3.0
crit=np.sqrt(2/3*gv*Lv*np.sin(bv)); print('critical v0',crit)
def rhs(t,y):
    X,U=y
    m=Lv-X
    # d/dt[(3/4)m U^2] + g sin b (L - X) U = 0, dm/dt=-U
    # (3/4)(-U)U^2 + (3/2) m U A + g sinb (L-X) U =0
    A=((3/4)*U**2 - gv*np.sin(bv)*(Lv-X))/(1.5*m)
    return [U,A]
ev=lambda t,y:y[1]; ev.terminal=True
so=solve_ivp(rhs,[0,100],[0,v0v],events=ev,rtol=1e-10,atol=1e-12)
xstop=so.y_events[0][0][0]
th=float(xs[0].subs({L:Lv,v0:v0v,g:gv,b:bv}))
print('numeric stop x',xstop,' formula',th, xs)
# horizontal time numeric
def rhs0(t,y):
    X,U=y; m=Lv-X; return [U,(0.75*U**2)/(1.5*m)]
ev2=lambda t,y:y[0]-Lv*0.999; ev2.terminal=True
so2=solve_ivp(rhs0,[0,100],[0,v0v],events=ev2,rtol=1e-11,atol=1e-12)
print('numeric t(0.999L)',so2.t_events[0][0],' formula',2*Lv/(3*v0v)*(1-0.001**1.5))
# numbers for problem: L=20, beta: sin b = 0.10, v0=2.8 -> check
vals={L:20,g:9.8,b:sp.asin(sp.Rational(1,10)),v0:sp.Rational(14,5)}
print('crit v0',float(sp.sqrt(sp.Rational(2,3)*9.8*20*0.1)))
print('stop x',[float(r.subs(vals)) for r in xs])
