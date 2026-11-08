# P31-1 ばね付きピストンと真空容器への自由膨張・かき混ぜ仕事
import sympy as sp
k,S,x0,n,R=sp.symbols('k S x0 n R',positive=True)
u,x1,x=sp.symbols('u x1 x',positive=True)   # u=C_V/R
a=2*x0   # V_B=2 S x0
T0=k*x0**2/(n*R)
def Tof(xx,open_):  # equilibrium temperature at piston height xx
    V=S*xx+(S*a if open_ else 0); P=k*xx/S
    return sp.simplify(P*V/(n*R))
def E(xx,open_):    # U + spring energy
    return n*u*R*Tof(xx,open_)+k*xx**2/2
# (1) free expansion: E conserved, x1=x0/2 -> u
usol=sp.solve(sp.Eq(E(x0,False),E(x0/2,True)),u)[0]; print('C_V/R=',usol)
T1=Tof(x0/2,True); print('T1/T0=',sp.simplify(T1/T0))
# also: general x1 eq for u=3/2 has root x0/2
print('roots x1:',sp.solve(sp.Eq(E(x0,False).subs(u,usol),E(x1,True).subs(u,usol)),x1))
# (2) stirring work to return to x0 with cock open
W=sp.simplify((E(x0,True)-E(x0/2,True)).subs(u,usol)); print('W=',W,'  T2/T0=',sp.simplify(Tof(x0,True)/T0))
# breakdown
dU=sp.simplify((n*u*R*(Tof(x0,True)-T1)).subs(u,usol)); dEs=sp.simplify(k*x0**2/2-k*(x0/2)**2/2); print('dU=',dU,' dEsp=',dEs)
# (3) heater with cock closed (quasi-static, P=kx/S), then open -> x0
Qp=sp.simplify((E(x0,True)-E(x0,False)).subs(u,usol)); print("Q'=",Qp)
xp=sp.solve(sp.Eq(E(x,False).subs(u,usol),E(x0,False).subs(u,usol)+Qp),x); print("x' =",xp)
# independent check of Q' by first law along quasi-static path: Q = dU + int P dV
xs=sp.symbols('xs',positive=True)
xpv=xp[0]
Qcheck=sp.simplify(n*usol*R*(Tof(xpv,False)-T0)+sp.integrate(k*xs/S*S,(xs,x0,xpv)))
print('Q via 1st law =',Qcheck)
