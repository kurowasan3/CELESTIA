# T20163-1: boat-borne "Michelson-Morley" for water waves
import sympy as sp, numpy as np
from scipy.optimize import brentq
d,V,u,f0,b,D=sp.symbols('d V u f0 beta D',positive=True)
tpar=d/(V-u)+d/(V+u); tperp=2*d/sp.sqrt(V**2-u**2)
print('t_par',sp.simplify(tpar))
# antiphase: f0*(tpar-tperp)=1/2 with D=d f0/V, beta=u/V
eq=2*D*(1/(1-b**2)-1/sp.sqrt(1-b**2))-sp.Rational(1,2)
s=sp.symbols('s',positive=True)
ss=sp.solve(sp.Eq(2*D*(s**2-s),sp.Rational(1,2)),s); print('s=',ss)
b2=sp.nsimplify(sp.expand(1-1/((2+sp.sqrt(5))/4)**2)); print('beta^2 (D=4)=',sp.simplify(b2),float(b2),' beta=',float(sp.sqrt(b2)))
print('check eq',float(eq.subs({D:4,b:sp.sqrt(b2)})))
# numeric geometry: water frame, boat at x=u t, plate1 line y=d (moving along itself), plate2 at x=u t + d
VV,uu,dd=1.0,0.33,1.0
# perpendicular arm: emitted at t=0 from origin, reflect at y=d at point (X,d); return to boat at time T: (uT,0)
def F(T): return 2*np.hypot(uu*T/2,dd)-VV*T
T1=brentq(F,0.1,100); print('num t_perp',T1,' th',2*dd/np.sqrt(VV**2-uu**2))
# emission angle from y axis: tan = (uT/2)/d
print('sin phi num',(uu*T1/2)/np.hypot(uu*T1/2,dd),' th',uu/VV)
# parallel arm
ta=brentq(lambda t: VV*t-(dd+uu*t),0,100); tb=brentq(lambda t: (dd+uu*ta)-VV*(t-ta)-uu*t,ta,100)
print('num t_par',tb,' th',2*dd*VV/(VV**2-uu**2))
# frequency: emission times t_e and t_e+dt arrive with same delay (translational invariance) -> f0. Check fixed front wall Doppler:
# wave emitted at t_e from x=u t_e toward wall at x=W (fixed): arrival back at boat
W=5.0
def arrive(te):
    xe=uu*te; tw=te+(W-xe)/VV
    return brentq(lambda t: (W-VV*(t-tw))-uu*t, tw, tw+100)
dt=1e-3; ratio=dt/(arrive(dt)-arrive(0)); print('fixed wall freq ratio num',ratio,' th',(VV+uu)/(VV-uu))
