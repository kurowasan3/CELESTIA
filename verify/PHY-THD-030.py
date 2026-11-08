# J061-1: boiling under a free piston: latent heat = Delta U + p Delta V
import sympy as sp
P,R,T0,M,L,p0,S,mw=sp.symbols('P R T_0 M L p_0 S m_w',positive=True)
# boiling: mass rate P/L ; volume rate (P/L)/M * R T0/p0
u1=P/L/M*R*T0/(p0*S); dh=mw/M*R*T0/(p0*S)
frac=sp.simplify(p0*S*u1/P)
# after: n=mw/M, Cp=4R (Cv=3R): dT/dt=P/(n*4R); dV/dt=nR dT/dt/p0
n=mw/M; dTdt=P/(n*4*R); u2=sp.simplify(n*R*dTdt/(p0*S))
print('u1=',u1,' dh=',dh,' frac=',frac,' u2=',u2,' u1/u2=',sp.simplify(u1/u2))
Lsol=sp.solve(sp.Eq(u1/u2, sp.Symbol('u1')/sp.Symbol('u2')),L); print('L=',Lsol)
vals={R:8.31,T0:373,M:1.8e-2,L:2.26e6}
print('frac=',float(frac.subs(vals)),' u1/u2=',float((u1/u2).subs(vals)))
# alternative: first law per mole: Q = L*M per mol; W = p0*v_gas = R T0; dU = L M - R T0
print('dU per kg =',float((L-R*T0/M).subs(vals)))
# independent check of u2 via first law: Q = dU + W ; dU = n*3R dT ; W = n R dT ; => dT = Q/(4nR); W/Q = 1/4
print('W/Q after =',sp.Rational(1,4), ' -> u2 = P/(4 p0 S) consistent:', sp.simplify(u2-P/(4*p0*S))==0)
vals2={R:8.3,T0:373,M:1.8e-2,L:2.3e6}
print('problem values: f=',float(frac.subs(vals2)),' u1/u2=',float((u1/u2).subs(vals2)))
