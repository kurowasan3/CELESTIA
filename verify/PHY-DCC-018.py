# S0510-1 verification: free-electron (Drude-type) model, two regimes and Joule heating
import sympy as sp, numpy as np
from scipy.optimize import brentq
m,e,E,l,u,n,kB,T=sp.symbols('m e E l u n k_B T',positive=True)
a=e*E/m
tau=l/u
vbar_low=a*tau/2
rho=sp.simplify(E/(n*e*vbar_low)); print('rho =',rho)
assert sp.simplify(rho-2*m*u/(n*e**2*l))==0
# Joule: energy lost per collision (1/2) m (a tau)^2, rate 1/tau, per volume n
P=n*sp.Rational(1,2)*m*(a*tau)**2/tau
assert sp.simplify(P-E**2/rho)==0; print('Joule power per volume = E^2/rho OK')
# temperature dependence
uT=sp.sqrt(3*kB*T/m); print('rho(2T)/rho(T)=',sp.simplify(rho.subs(u,uT).subs(T,2*T)/rho.subs(u,uT)))
# high field
vbar_high=sp.sqrt(2*a*l)/2; print('vbar_high=',sp.simplify(vbar_high))
Ex=sp.solve(sp.Eq(vbar_low,vbar_high),E); print('crossover E =',Ex)
Exs=sp.simplify(Ex[0].subs(u,uT)); print('in terms of T:',Exs)
assert sp.simplify(Exs-6*kB*T/(e*l))==0
num=float(Exs.subs({kB:1.38e-23,T:300,e:1.6e-19,l:4.0e-8})); print('E_x numeric = %.3g V/m'%num)
# Direction-averaged numerical model (dimensionless m=e=l=u=1): after each collision speed u in random 3D direction,
# flight ends when path length reaches l. Estimate drift = <displacement along field>/<time>.
rng=np.random.default_rng(1)
def flight(Ev,N=20000):
    cz=(np.arange(N)+0.5)/N*2-1; s=np.sqrt(1-cz**2)  # deterministic midpoint quadrature over isotropic directions
    vz0=cz; vt=s   # u=1
    disp=0.0; tt=0.0
    # path length L(t)=int sqrt(vt^2+(vz0+Ev t)^2) dt ; solve L(t)=1 numerically per sample (vectorized bisection)
    lo=np.zeros(N); hi=np.full(N,10.0)
    def L(t):
        # analytic arc length of hyperbola-like: integral sqrt(A^2 + w^2) dw / Ev, w = vz0+Ev t
        w1=vz0+Ev*t; w0=vz0; A=vt
        F=lambda w: 0.5*(w*np.sqrt(A**2+w**2)+A**2*np.arcsinh(np.where(A>0,w/np.maximum(A,1e-300),0)))
        return (F(w1)-F(w0))/Ev
    for _ in range(80):
        mid=(lo+hi)/2; f=L(mid)-1.0
        hi=np.where(f>0,mid,hi); lo=np.where(f>0,lo,mid)
    t=(lo+hi)/2
    z=vz0*t+0.5*Ev*t**2
    return z.mean()/t.mean()
for Ev in [1e-3,1e-2,1e2,1e3]:
    vd=flight(Ev); low=Ev/2; high=np.sqrt(Ev/2)
    print(f'E={Ev:g}: MC drift={vd:.5g}, low-field formula={low:.5g}, high-field formula={high:.5g}')
# constant-path-length model: low field drift is ohmic (proportional to E) but with factor 2/3 relative to the
# constant-free-time model used in the problem (flight time correlates with initial direction); high-field limit matches exactly.
r1=flight(1e-3)/1e-3; r2=flight(1e-2)/1e-2
assert abs(r1-r2)<1e-3 and abs(r1/(0.5)-2/3)<0.01
assert abs(flight(1e3)/np.sqrt(1e3/2)-1)<0.01
# crossover check: constant-path model drift at E=Ex(=2 in these units) lies between the two regimes
print('E=2 (crossover):',flight(2.0),' formulas:',1.0,1.0)
print('ALL OK')
