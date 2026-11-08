# P51-1: grating with two slits per period (amplitudes 1 and r, offset s)
import sympy as sp, numpy as np
r,c,m,u=sp.symbols('r c m u',real=True)
Im=1+r**2+2*r*sp.cos(2*sp.pi*m*u)   # u=s/d
print('I0=',sp.factor(Im.subs(m,0)))
# s=d/3, r=1
print('s=d/3,r=1: I_m/I_0 for m=0..6:',[sp.nsimplify(Im.subs({r:1,u:sp.Rational(1,3),m:k})/4) for k in range(7)])
# inverse: I1/I0=5/9, I2/I0=1/9
eqs=[sp.Eq(1+r**2+2*r*c,sp.Rational(5,9)*(1+r)**2), sp.Eq(1+r**2+2*r*(2*c**2-1),sp.Rational(1,9)*(1+r)**2)]
sols=sp.solve(eqs,[r,c],dict=True); print('solutions (r,c):',sols)
# numeric grating with N periods, exact sum at order angles
d=2e-6; lam=0.5e-6; N=200
for (rr,uu) in [(0.5,0.25),(0.5,0.75)]:
    ys=np.concatenate([np.arange(N)*d, np.arange(N)*d+uu*d]); amps=np.concatenate([np.ones(N),rr*np.ones(N)])
    def I(sinth): return abs(np.sum(amps*np.exp(-2j*np.pi*ys*sinth/lam)))**2
    I0=I(0); print('u=',uu,[round(I(k*lam/d)/I0,4) for k in range(4)])
    # check maxima positions: scan near m=1
    th=np.linspace(0.24,0.26,20001); vals=[I(t) for t in th]; print('  peak sin at',th[int(np.argmax(vals))],'expected',lam/d)
