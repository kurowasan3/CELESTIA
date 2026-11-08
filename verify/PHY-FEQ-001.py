# P01-1 capstan on cylinder at top of incline: wrap angle, equilibrium range, kinetic phase
import sympy as sp, numpy as np
M,m,g,al,mu,mu1,mup,mu1p,d=sp.symbols('M m g alpha mu mu1 muprime mu1prime d',positive=True)
phi=sp.pi/2+al
# geometry check: angle between direction up-slope (cos a, sin a) and down (0,-1)
a0=0.4; u=np.array([np.cos(a0),np.sin(a0)]); w=np.array([0,-1.])
print('turn angle',np.arccos(u@w), np.pi/2+a0)
# discrete check of capstan factor: n pegs each turning phi/n, at each peg T_hi <= T_lo*(1+mu*tan(dphi/2))/(1-mu*tan(dphi/2))
def factor(muv,ph,n):
    t=np.tan(ph/n/2); return ((1+muv*t)/(1-muv*t))**n
for n in [3,10,100,10000]: print(n,factor(0.3,np.pi/2+a0,n))
print('exp',np.exp(0.3*(np.pi/2+a0)))
mmax=M*(sp.sin(al)+mu1*sp.cos(al))*sp.exp(mu*phi)
mmin=M*(sp.sin(al)-mu1*sp.cos(al))*sp.exp(-mu*phi)
# kinetic phase: m moves down, block moves up slope
a,T1,T2=sp.symbols('a T1 T2')
sol=sp.solve([sp.Eq(m*a,m*g-T1),sp.Eq(M*a,T2-M*g*(sp.sin(al)+mu1p*sp.cos(al))),sp.Eq(T1,T2*sp.exp(mup*phi))],[a,T1,T2],dict=True)[0]
k=sp.exp(-mup*phi)
aform=g*(m*k-M*(sp.sin(al)+mu1p*sp.cos(al)))/(M+m*k)
print('a check',sp.simplify(sol[a]-aform))
print('T1-T2',sp.simplify(sol[T1]-sol[T2]))
# energy check: over distance d: m g d - M g d sin a = 1/2(M+m)v^2 + heat_floor + heat_cyl
v2=2*sol[a]*d
heat_cyl=(sol[T1]-sol[T2])*d; heat_incl=mu1p*M*g*sp.cos(al)*d
print('energy',sp.simplify(m*g*d-M*g*d*sp.sin(al)-sp.Rational(1,2)*(M+m)*v2-heat_cyl-heat_incl))
# numeric example for answers
vals={M:1,g:9.8,al:sp.pi/6,mu:0.3,mu1:0.2}
print('mmax,mmin',sp.N(mmax.subs(vals)),sp.N(mmin.subs(vals)))
# positivity: T1>0 automatically; also need a>0 for m > M(s+mu1'c)e^{mu' phi}; since mu'<mu and mu1'<mu1, mmax > threshold -> a>0 just above mmax
