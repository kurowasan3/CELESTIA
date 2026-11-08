# S0313-1: capillary-gravity waves via Bernoulli in comoving frame
import sympy as sp, numpy as np
from scipy.optimize import brentq
rho,g,s,a,k,v,w,U=sp.symbols('rho g sigma a k v omega U',positive=True)
# (1) Bernoulli along surface streamline in wave frame: crest p0+s a k^2, trough p0 - s a k^2
eq=sp.Eq(s*a*k**2+sp.Rational(1,2)*rho*(v-a*w)**2+rho*g*a, -s*a*k**2+sp.Rational(1,2)*rho*(v+a*w)**2-rho*g*a)
sol=sp.solve(eq.subs(w,v*k),v)
print('(1) v =',sol)
v2=g/k+s*k/rho
# (2) minimum
km=sp.solve(sp.diff(v2,k),k)[0]; vmin=sp.sqrt(v2.subs(k,km))
print('(2) k_m=',km,' vmin=',sp.simplify(vmin))
num={rho:1.0e3,g:9.8,s:7.3e-2}
print('    lambda_m=%.4g m vmin=%.4g m/s'%(float(2*sp.pi/km.subs(num)),float(vmin.subs(num))))
# (3) group velocity
om=sp.sqrt(g*k+s*k**3/rho)
vg=sp.diff(om,k); vp=om/k
print('(3) vg/v at km:',sp.simplify((vg/vp).subs(k,km)))
print('    limits:',sp.limit(vg/vp,k,0),sp.limit(vg/vp,k,sp.oo))
print('    vg/v general =',sp.simplify(vg/vp))
# (4) U=0.30 m/s stream: solve v(k)=U
Un=0.30; gn,sn,rn=9.8,7.3e-2,1.0e3
f=lambda kk: gn/kk+sn*kk/rn-Un**2
kmn=float(km.subs(num))
k1=brentq(f,1e-3,kmn); k2=brentq(f,kmn,1e6)
for kk in (k1,k2):
    vgn=(gn+3*sn*kk**2/rn)/(2*np.sqrt(gn*kk+sn*kk**3/rn))
    print('    k=%.4g lambda=%.4g cm  vg=%.4g  (U=%.2f) -> %s'%(kk,200*np.pi/kk,vgn,Un,'upstream' if vgn>Un else 'downstream'))
# quadratic form: s k^2/rho - U^2 k + g = 0
kq=sp.solve(sp.Eq(s*k**2/rho-U**2*k+g,0),k)
print('    quadratic roots:',[float(r.subs(num).subs(U,Un)) for r in kq])
# numerical Bernoulli check of (1) independent: energy method with given numbers
