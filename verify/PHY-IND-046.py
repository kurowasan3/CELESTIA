# S0525-1 verification: RL ring in linearly ramped external flux
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
t,T,L,R,P0=sp.symbols('t T L R Phi0',positive=True)
tau=L/R; a=P0/T
I=sp.Function('I')
# circuit: R I = -dPhi_e/dt - L dI/dt  (0<t<T)
sol=sp.dsolve(sp.Eq(R*I(t),-a-L*I(t).diff(t)),I(t),ics={I(0):0}).rhs
Iexp=-(P0/(R*T))*(1-sp.exp(-t/tau))
print('(1)',sp.simplify(sol-Iexp))
IT=Iexp.subs(t,T)
PhiT=P0+L*IT
print('(2) deficit',sp.simplify(P0-PhiT-(P0*tau/T)*(1-sp.exp(-T/tau))))
# heat: direct integral during ramp + L IT^2/2 after
H1=sp.integrate(R*Iexp**2,(t,0,T))
H=sp.simplify(H1+L*IT**2/2)
Hexp=P0**2/(R*T**2)*(T-tau*(1-sp.exp(-T/tau)))
print('(4)',sp.simplify(H-Hexp))
# energy-method: work of emf (-a) times ramp charge
q=sp.integrate(Iexp,(t,0,T))
print('(4) emf-work method',sp.simplify(-a*q-Hexp))
print('(5) T->0',sp.limit(Hexp,T,0),' T->oo leading',sp.limit(Hexp*T,T,sp.oo))
# numeric check incl. charge after T
Lv,Rv,P0v,Tv=0.3,2.0,1.5,0.4
def f(tt,y):
    dPe=P0v/Tv if tt<Tv else 0.0
    I_=y[0]
    return [(-dPe-Rv*I_)/Lv, I_, Rv*I_**2]
s=solve_ivp(f,[0,20],[0,0,0],rtol=1e-10,atol=1e-12,max_step=0.01)
print('num Q',s.y[1,-1],' expect',-P0v/Rv)
Hn=float(Hexp.subs({L:Lv,R:Rv,P0:P0v,T:Tv}))
print('num H',s.y[2,-1],' expect',Hn)
# monotonic decrease of H in T
Tg=np.linspace(0.01,10,500); f_=sp.lambdify(T,Hexp.subs({L:Lv,R:Rv,P0:P0v}))
print('monotone decreasing:',np.all(np.diff(f_(Tg))<0))
