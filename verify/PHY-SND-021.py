# S0305-1 verification: ring-tube acoustic interferometer with circulating air / rotating tube
import sympy as sp, numpy as np
from scipy.optimize import brentq
f,l,V,u,u1,u2=sp.symbols('f ell V u u1 u2',positive=True)
lam_cw=(V+u)/f; lam_ccw=(V-u)/f
dphi=sp.simplify(2*sp.pi*(l/lam_ccw-l/lam_cw)); print('dphi =',sp.factor(dphi))
dphi_t=sp.simplify(2*sp.pi*f*(l/(V-u)-l/(V+u))); print('time method equal:',sp.simplify(dphi-dphi_t)==0)
# first silence
us=sp.solve(sp.Eq(dphi,sp.pi),u); print('u_first =',us)
# inverse
Vs,F=sp.symbols('Vs F',positive=True)
sol=sp.solve([sp.Eq(4*F*u1,Vs**2-u1**2),sp.Eq(4*F*u2,3*(Vs**2-u2**2))],[Vs,F],dict=True)
print('inverse:',[ {k:sp.factor(v) for k,v in s.items()} for s in sol])
# numeric round trip
fn,ln,Vn=500.0,2.0,340.0
g=lambda uu,k: 4*fn*ln*uu/(Vn**2-uu**2)-(2*k+1)
U1=brentq(lambda x:g(x,0),1e-9,Vn-1e-6); U2=brentq(lambda x:g(x,1),1e-9,Vn-1e-6)
Vrec=np.sqrt(U1*U2*(3*U2-U1)/(3*U1-U2)); print('U1,U2',U1,U2,' V recovered',Vrec,' f*l recovered',(Vrec**2-U1**2)/(4*U1))
# (3) lab frame: air at rest, tube (source P, mic Q) rotates; rim speed un.
# arc coordinate s along ring (circumference 2l) in lab; P at s=un*t, Q at s=l+un*t.
un=30.0; t_obs=1.0
# CW-going sound (direction +s): emitted at te from P(te), arrives Q at t_obs: V*(t_obs-te) = (l+un*t_obs) - un*te  (mod 2l, take first)
te_cw=brentq(lambda te: Vn*(t_obs-te)-((ln+un*t_obs)-un*te),t_obs-1,t_obs)
# CCW (-s): V*(t_obs-te) = un*te + 2l - (l+un*t_obs)... distance going backwards from P(te) to Q(t_obs): P - Q + 2l
te_ccw=brentq(lambda te: Vn*(t_obs-te)-((un*te)-(ln+un*t_obs)+2*ln),t_obs-1,t_obs)
dphi_lab=2*np.pi*fn*(te_ccw-te_cw)  # source phase 2 pi f te (source clock = lab clock, Galilean)
print('lab-frame dphi',dphi_lab,' tube-frame formula',4*np.pi*fn*ln*un/(Vn**2-un**2))
# received frequency: d(te)/dt_obs for each path should be 1
h=1e-3
te_cw2=brentq(lambda te: Vn*(t_obs+h-te)-((ln+un*(t_obs+h))-un*te),t_obs+h-1,t_obs+h)
print('received freq factor CW',(te_cw2-te_cw)/h)
