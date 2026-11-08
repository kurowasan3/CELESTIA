# T20061-1 pulsar planet: timing vs Doppler(pulse interval)
import sympy as sp
import numpy as np

G, M, m, T, c, tau0, eps = sp.symbols('G M m T c tau0 epsilon', positive=True)
d = (G*(M+m)*T**2/(4*sp.pi**2))**sp.Rational(1, 3)
# check Kepler: two-body circular orbit, relative motion with reduced mass
om = 2*sp.pi/T
assert sp.simplify(om**2*d - G*(M+m)/d**2) == 0
a = m*d/(M+m)
V = sp.simplify(2*sp.pi*a/T)
V_claim = m*(2*sp.pi*G/T)**sp.Rational(1, 3)*(M+m)**sp.Rational(-2, 3)
print('V check', sp.simplify(V/V_claim))
a_claim = m*(G*T**2/(4*sp.pi**2))**sp.Rational(1, 3)*(M+m)**sp.Rational(-2, 3)
print('a check', sp.simplify(a/a_claim))
# also check the pulsar's own equation of motion: M*om^2*a = G M m / d^2
print('EOM pulsar', sp.simplify(M*om**2*a - G*M*m/d**2))

# minimal masses (m<<M): timing a/c=tau0 ; Doppler V/c=eps
mt = tau0*c*M**sp.Rational(2, 3)*(4*sp.pi**2/(G*T**2))**sp.Rational(1, 3)
mD = eps*c*M**sp.Rational(2, 3)*(T/(2*sp.pi*G))**sp.Rational(1, 3)
a0 = a_claim.subs(m+M, M).subs(M+m, M)
a0 = m*(G*T**2/(4*sp.pi**2))**sp.Rational(1, 3)*M**sp.Rational(-2, 3)
V0 = m*(2*sp.pi*G/T)**sp.Rational(1, 3)*M**sp.Rational(-2, 3)
print('mt check', sp.simplify(sp.solve(sp.Eq(a0/c, tau0), m)[0]/mt))
print('mD check', sp.simplify(sp.solve(sp.Eq(V0/c, eps), m)[0]/mD))
print('ratio mt/mD', sp.simplify(mt/mD))

# numerics
vals = {G: 6.7e-11, M: 2.8e30, c: 3.0e8, tau0: 1.0e-6, eps: 1.0e-11}
Tc = float((2*sp.pi*tau0/eps).subs(vals))
print('Tc [s]', Tc, ' days', Tc/86400)
Ty = 3.2e7
mtn = float(mt.subs(vals).subs(T, Ty)); mDn = float(mD.subs(vals).subs(T, Ty))
print('T=3.2e7: m_timing', mtn, ' /Earth(6.0e24)', mtn/6.0e24, ' m_Doppler', mDn, mDn/6.0e24)

# numerical simulation of pulse arrival: pulses emitted with P0, star line-of-sight x=a sin(wt)
mm = 6.0e24; MM = 2.8e30; Gv = 6.7e-11; cv = 3.0e8; Tv = 3.2e7
dd = (Gv*(MM+mm)*Tv**2/(4*np.pi**2))**(1/3); aa = mm*dd/(MM+mm); w = 2*np.pi/Tv
P0 = 1.0e3  # coarse emission interval for test
te = np.arange(0, 2*Tv, P0)
tarr = te + aa*np.sin(w*te)/cv
dt = np.diff(tarr)
print('timing amplitude a/c', aa/cv, ' sim', (tarr-te).max())
print('interval frac max', (dt.max()-P0)/P0, ' V/c', w*aa/cv)
