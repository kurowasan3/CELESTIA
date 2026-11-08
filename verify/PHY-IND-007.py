# N42-1: Faraday disk as motor with battery then generator charging capacitor
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
B,a,J,R,E,C,t=sp.symbols('B a J R E C t',positive=True)
K=B*a**2/2
# (1) emf by integrating Lorentz-force-balanced field E=w B r
w,rr=sp.symbols('omega r',positive=True)
print('emf =',sp.integrate(w*B*rr,(rr,0,a)))
# (2) motor: J dw/dt = K (E - K w)/R
wf=sp.Function('w')
sol=sp.dsolve(sp.Eq(J*wf(t).diff(t),K*(E-K*wf(t))/R),wf(t),ics={wf(0):0})
print('w(t) =',sp.simplify(sol.rhs))
winf=E/K; print('w_inf =',sp.simplify(winf),' tau =',sp.simplify(J*R/K**2))
Ws=sp.integrate(E*K*(E-K*sol.rhs)/R/K*K/K*1,(t,0,sp.oo))  # battery power E*I, I=(E-Kw)/R
Ws=sp.integrate(E*(E-K*sol.rhs)/R,(t,0,sp.oo)); Heat=sp.integrate(R*((E-K*sol.rhs)/R)**2,(t,0,sp.oo))
print('battery work',sp.simplify(Ws),' heat',sp.simplify(Heat),' ratio',sp.simplify(Heat/Ws))
# (3) capacitor stage: conservation J w + K q = J w_inf ; final K w = q/C
wfin=sp.symbols('w_f'); wfs=sp.solve(sp.Eq(J*wfin+K*C*K*wfin,J*winf),wfin)[0]
print('w_f =',sp.simplify(wfs)); qf=C*K*wfs; print('q_f =',sp.simplify(qf))
H2=sp.simplify(J*winf**2/2-J*wfs**2/2-qf**2/(2*C)); print('heat2 =',sp.factor(H2))
Ceq=J/K**2; print('heat2 vs CeqC/(Ceq+C) E^2/2:',sp.simplify(H2-Ceq*C/(Ceq+C)*E**2/2))
# numeric ODE
vals={B:0.5,a:0.2,J:2e-3,R:0.05,E:1.5,C:0.02}
Kv=float(K.subs(vals)); Jv=vals[J];Rv=vals[R];Ev=vals[E];Cv=vals[C]
s1=solve_ivp(lambda t,y:[Kv*(Ev-Kv*y[0])/Rv/Jv,Rv*((Ev-Kv*y[0])/Rv)**2],[0,50*Jv*Rv/Kv**2],[0,0],rtol=1e-10,atol=1e-12)
print('num w_end',s1.y[0,-1],' formula',float(winf.subs(vals)),' heat1',s1.y[1,-1],' 1/2Jw^2',0.5*Jv*s1.y[0,-1]**2)
def rhs(t,y):
    ww,q,h=y; I=(Kv*ww-q/Cv)/Rv
    return [-Kv*I/Jv,I,Rv*I**2]
s2=solve_ivp(rhs,[0,5],[s1.y[0,-1],0,0],rtol=1e-11,atol=1e-13)
print('num w_f',s2.y[0,-1],' formula',float(wfs.subs(vals)),' q',s2.y[1,-1],float(qf.subs(vals)),' heat2',s2.y[2,-1],float(H2.subs(vals)))
