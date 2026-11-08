# O20233-1 balloon with two gases (mono A, di B), isothermal atmosphere T0
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
n,M,m,T0,p0,R,g,a=sp.symbols('n M m T0 p0 R g a',positive=True)
# buoyancy of bag at ambient pressure p, gas temps TA,TB, air T0: rho_air*g*V
def buoy(TA,TB,p): return p*M/(R*T0)*g*(n*R*(TA+TB)/p)
T1=sp.solve(sp.Eq(buoy(sp.Symbol('T'),sp.Symbol('T'),p0),m*g),sp.Symbol('T'))[0]
print('T1=',sp.simplify(T1))
Q1=sp.simplify(n*(sp.Rational(5,2)+sp.Rational(7,2))*R*(T1-T0)); print('Q1=',sp.factor(Q1))
T2=sp.Symbol('T2',positive=True)
TA=T2*a**sp.Rational(2,5); TB=T2*a**sp.Rational(2,7)   # pV^g const -> T p^{(1-g)/g}
# check adiabatic exponent: gamma=5/3 -> (g-1)/g=2/5 ; 7/5 -> 2/7
T2sol=sp.solve(sp.Eq(buoy(TA,TB,a*p0),m*g),T2)[0]; print('T2=',sp.simplify(T2sol))
TAv=TA.subs(T2,T2sol); TBv=TB.subs(T2,T2sol)
Tf=(5*TAv+7*TBv)/12
S=sp.simplify(buoy(Tf,Tf,a*p0)-m*g); 
S_ans=m*g*(a**sp.Rational(2,7)-a**sp.Rational(2,5))/(6*(a**sp.Rational(2,5)+a**sp.Rational(2,7)))
print('tension check:',sp.simplify(S-S_ans))
# release: combined adiabatic gamma=3/2 => T ∝ p^{1/3}; float when 2T=2T1
pp=a*p0*(T1/Tf)**3
dh=R*T0/(M*g)*sp.log(a*p0/pp)
dh_ans=3*R*T0/(M*g)*sp.log((5*a**sp.Rational(2,5)+7*a**sp.Rational(2,7))/(6*(a**sp.Rational(2,5)+a**sp.Rational(2,7))))
print('dh check:',sp.simplify(sp.expand_log(dh-dh_ans,force=True)))
# numeric: integrate first law with common T during slow ascent, find float pressure
vals={n:1.0,M:0.029,m:0.09,T0:280.0,p0:1.0e5,R:8.31,g:9.8,a:0.6}
T1n=float(T1.subs(vals)); Tfn=float(Tf.subs(vals)); print('T1,Tf,TA,TB',T1n,Tfn,float(TAv.subs(vals)),float(TBv.subs(vals)))
# state (T) as function of p: dU=-p dV ; U=4nRT (3/2+5/2), V=2nRT/p
def rhs(p,y):
    T=y[0]; # 4R dT = -p d(2RT/p) = -2R dT + 2R T dp/p
    return [2*T/p/6]
pa=0.6e5
sol=solve_ivp(rhs,[pa,0.2e5],[Tfn],dense_output=True,rtol=1e-11,atol=1e-9)
f=lambda p: sol.sol(p)[0]-T1n
pnum=brentq(f,0.2e5,pa); print('p\' numeric',pnum,'formula',float(pp.subs(vals)))
print('dh',float(dh_ans.subs(vals)), 'numeric', 8.31*280/(0.029*9.8)*np.log(pa/pnum))
# direct check of tension from volumes
V_before=8.31*(float(TAv.subs(vals))+float(TBv.subs(vals)))/pa; V_after=8.31*2*Tfn/pa
print('tension numeric', pa*0.029/(8.31*280)*9.8*(V_after-V_before), float(S_ans.subs(vals)))
print('T_B>T_A ?', float(TBv.subs(vals))>float(TAv.subs(vals)), ' T1>T0?',T1n>280)
