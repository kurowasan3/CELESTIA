# N53-1 180-degree magnetic focusing: landing distance 2R cos(phi); resolution
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
M,Q,E,B1,B2,phi,dl,eps=sp.symbols('M Q E B_1 B_2 phi delta epsilon',positive=True)
v=E/B1; R=M*v/(Q*B2); D0=2*R
# circle: enter at origin with velocity (sin phi, cos phi) into y<0? set boundary y=0, field region y<0, ion enters moving -y
# geometric: center at distance R perpendicular to velocity; chord on y=0 line = 2R cos(phi)
print('D0=',sp.simplify(D0))
# resolution conditions
c1=sp.acos(sp.Rational(35,37)); print('delta max (no speed spread) deg', float(c1*180/sp.pi))
print('approx sqrt(2*2/35) deg', float(sp.sqrt(sp.Rational(4,35))*180/sp.pi))
# numeric integration of landing point vs phi
Mn,Qn,Bn,vn=6e-26,1.6e-19,0.5,1e5
Rn=Mn*vn/(Qn*Bn)
def land(ph):
    def f(t,s):
        vel=np.r_[s[2],s[3],0]; F=Qn*np.cross(vel,[0,0,-Bn]); return np.r_[vel[:2],F[:2]/Mn]
    ev=lambda t,s: s[1]; ev.terminal=True; ev.direction=1
    sol=solve_ivp(f,[0,1e-3],[0,0,vn*np.sin(ph),-vn*np.cos(ph)],events=ev,rtol=1e-11,atol=1e-14,first_step=1e-12)
    return sol.y_events[0][0][0]
for ph in [0,0.2,-0.2]:
    print(ph, land(ph), 2*Rn*np.cos(ph))

# (5) new: general condition with relative speed spread eps and mass ratio mu=dM/M:
# (1+mu)(1-eps)cos(delta) > 1+eps  ->  first order: delta^2 < 2(mu-2 eps)
mu,lam=sp.symbols('mu lambda',positive=True)
expr=(1+mu)*(1-eps)*sp.cos(dl)-(1+eps)
# scale all small quantities: mu->lam*mu, eps->lam*eps, delta->sqrt(lam)*delta, expand to first order in lam
first=sp.series(expr.subs({mu:lam*mu,eps:lam*eps,dl:sp.sqrt(lam)*dl}),lam,0,2).removeO()
print('first order:',sp.expand(first/lam))
assert sp.simplify(first/lam-(mu-2*eps-dl**2/2))==0
# exact vs approximate threshold for some small values
for muv,ev in [(2/35,0.005),(0.01,0.002),(0.003,0.001)]:
    ex=np.arccos((1+ev)/((1+muv)*(1-ev))); ap=np.sqrt(2*(muv-2*ev))
    print(f'mu={muv:.4f} eps={ev}: exact delta_max={ex:.4f}  approx={ap:.4f}')
# no delta works when eps >= mu/2 (exact: (1+mu)(1-eps) <= 1+eps  <=> eps >= mu/(2+mu) ~ mu/2)
print('exact threshold eps* =',sp.solve(sp.Eq((1+mu)*(1-eps),1+eps),eps),' ~ mu/2 to first order')
# order of spreading: landing distance vs phi (2nd order) and vs speed (1st order)
print('dD/dphi at 0:',sp.diff(D0*sp.cos(phi),phi).subs(phi,0),'  dD/dv relative:',sp.simplify(sp.diff(D0,E)*E/D0))
for ph in [0.05,0.1]:
    print('numeric landing shift for phi=',ph,':',(2*Rn-land(ph))/(2*Rn),' vs phi^2/2=',ph**2/2)
