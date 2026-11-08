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
# (4) 新設問：棒に対して静止する2つの波の波長の積 = lambda_m^2，上流/下流の選択，U->vmin で合流
kq=sp.solve(sp.Eq(s*k**2/rho-U**2*k+g,0),k)
prod=sp.simplify(kq[0]*kq[1]); print('(4a) k1*k2 =',prod,' k_m^2 =',sp.simplify(km**2))
assert sp.simplify(prod-km**2)==0
gn,sn,rn=9.8,7.3e-2,1.0e3
kmn=float(km.subs(num)); lam_m=2*np.pi/kmn; vminn=float(vmin.subs(num))
for Un in [0.24,0.30,0.5,1.0]:
    f=lambda kk: gn/kk+sn*kk/rn-Un**2
    k1=brentq(f,1e-6,kmn); k2=brentq(f,kmn,1e7)
    l1,l2=2*np.pi/k1,2*np.pi/k2
    vg=lambda kk:(gn+3*sn*kk**2/rn)/(2*np.sqrt(gn*kk+sn*kk**3/rn))
    print(f'U={Un}: lam1={l1*100:.4g} cm lam2={l2*100:.4g} cm  lam1*lam2/lam_m^2={l1*l2/lam_m**2:.12f}  vg(k1)-U={vg(k1)-Un:+.4f} vg(k2)-U={vg(k2)-Un:+.4f}')
    assert abs(l1*l2/lam_m**2-1)<1e-9 and vg(k2)>Un and vg(k1)<Un
# U -> vmin+ で両者が lam_m に合流し，vg - U -> 0
for eps_ in [1e-2,1e-4,1e-6]:
    Un=vminn*(1+eps_)
    f=lambda kk: gn/kk+sn*kk/rn-Un**2
    k1=brentq(f,1e-6,kmn); k2=brentq(f,kmn,1e7)
    vg=lambda kk:(gn+3*sn*kk**2/rn)/(2*np.sqrt(gn*kk+sn*kk**3/rn))
    print(f'U=vmin(1+{eps_}): lam1/lam_m={2*np.pi/k1/lam_m:.5f} lam2/lam_m={2*np.pi/k2/lam_m:.5f}  |vg-U| = {abs(vg(k1)-Un):.2e}, {abs(vg(k2)-Un):.2e}')
# 図の例 U=0.30 m/s
Un=0.30; f=lambda kk: gn/kk+sn*kk/rn-Un**2
k1=brentq(f,1e-3,kmn); k2=brentq(f,kmn,1e6)
print('figure example U=0.30: lambda1=%.3g cm lambda2=%.3g cm lambda_m=%.3g cm'%(200*np.pi/k1,200*np.pi/k2,100*lam_m))
print('PHY-WAV-019 (4) OK')
