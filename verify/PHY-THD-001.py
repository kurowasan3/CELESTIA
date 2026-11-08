# M47-1: open vs sealed hot-air balloon in isothermal atmosphere
import sympy as sp, numpy as np
from scipy.optimize import brentq
V,W,T0,T,rho0,P0,g,dw=sp.symbols('V W T0 T rho0 P0 g Delta_w',positive=True)
rhoh=sp.symbols('rho_h',positive=True)
T1=sp.solve(sp.Eq(rho0*V*(1-T0/T),W),T)[0]; print('T1',sp.simplify(T1))
rh=sp.solve(sp.Eq(rhoh*V*(1-T0/T),W),rhoh)[0]; print('rho_h',sp.simplify(rh))
# sealed: inner mass
min_=rh*T0/T*V
# sealed after cooling: weight unchanged -> balance unchanged
print('sealed balance residual',sp.simplify(rh*V-(W+min_)))
# small-change approx: drho = -(rho0 g/P0) rho dh
k=rho0*g/P0
dh_closed=sp.simplify((dw/V)/(k*rh)); dh_open=sp.simplify((rh*dw/W)/(k*rh))
print('dh_closed',dh_closed,' dh_open',dh_open,' ratio',sp.simplify(dh_open/dh_closed))
# exact numeric check with exponential atmosphere
vals=dict(V=2000.,W=600.,T0=280.,T=380.,rho0=1.20,P0=1.0e5,g=9.8,dw=5.0)
Vn,Wn,T0n,Tn,r0,P0n,gn,dwn=[vals[x] for x in ['V','W','T0','T','rho0','P0','g','dw']]
kk=r0*gn/P0n
rho=lambda h: r0*np.exp(-kk*h)
# open equilibrium
h_open0=brentq(lambda h: rho(h)*Vn*(1-T0n/Tn)-Wn,0,2e4)
h_open1=brentq(lambda h: rho(h)*Vn*(1-T0n/Tn)-(Wn-dwn),0,2e4)
rh_n=rho(h_open0); mi=rh_n*T0n/Tn*Vn
h_cl1=brentq(lambda h: rho(h)*Vn-(Wn-dwn+mi),0,2e4)
print('h0',h_open0,'open rise',h_open1-h_open0,'closed rise',h_cl1-h_open0)
sub={V:Vn,W:Wn,T0:T0n,T:Tn,rho0:r0,P0:P0n,g:gn,dw:dwn}
print('formula open',float(dh_open.subs(sub)),'closed',float(dh_closed.subs(sub)))
print('ratio exact',(h_open1-h_open0)/(h_cl1-h_open0),'formula',Tn/(Tn-T0n))
