# T20002-1 verification: parallel-strip line, electric attraction vs magnetic repulsion
import sympy as sp
e0, mu0, V, R, h, w, l = sp.symbols('epsilon0 mu0 V R h w ell', positive=True)
c = 1/sp.sqrt(e0*mu0)
sigma = e0*V/h                       # field between plates V/h = sigma/e0
Q = sigma*w*l
FE = Q*sigma/(2*e0)                  # force on upper plate from lower plate field sigma/(2e0)
I = V/R
K = I/w
FB = K*(mu0*K/2)*w*l                 # force per area K * B_other
print('F_E =', sp.simplify(FE))
print('F_B =', sp.simplify(FB))
R0 = sp.solve(sp.Eq(FE, FB), R)
print('R0 =', R0, ' == mu0 c h/w ?', sp.simplify(R0[0] - mu0*c*h/w) == 0)
UE = sp.Rational(1,2)*(e0*w*l/h)*V**2
Bin = mu0*K
UB = Bin**2/(2*mu0)*w*h*l
print('U_B/U_E =', sp.simplify(UB/UE), '  F_B/F_E =', sp.simplify(FB/FE))
print('ratio equals (R0/R)^2 :', sp.simplify(UB/UE - (R0[0]/R)**2) == 0)
# virtual work check: F_E = dU_E/dh at const Q ; F_B: at const I magnetic force = +dU_B/dh
hh = sp.symbols('hh', positive=True)
UE_Q = (Q**2/(2*(e0*w*l/hh)))
print('virtual work F_E:', sp.simplify(sp.diff(UE_Q, hh).subs(hh, h) - FE) == 0)
UB_I = mu0*I**2*hh*l/(2*w)
print('virtual work F_B:', sp.simplify(sp.diff(UB_I, hh).subs(hh, h) - FB) == 0)
# numerics
mu0v = 4*sp.pi*1e-7
print('mu0 c (c=3.0e8) =', sp.N(mu0v*3.0e8))
print('R0 for h/w=0.020:', sp.N(mu0v*3.0e8*0.020))
print('c from R0=6.0, h/w=0.016:', sp.N(6.0/(mu0v*0.016)))
print('check (4) numbers: V=100V, w=5cm, h=1mm, l=1m, R=R0 -> F_E =', sp.N(FE.subs({e0:8.85e-12, V:100, h:1e-3, w:0.05, l:1})))
