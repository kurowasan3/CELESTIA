# K20213-1 verification: gravity-induced phase in neutron interferometer with lambda/2 contamination
import sympy as sp, numpy as np
from scipy.optimize import minimize_scalar
m,g,h,l,s,lam,al,v0=sp.symbols('m g h ell s lambda_0 alpha v_0',positive=True)
# exact wavelength on DC side
v1=sp.sqrt(v0**2-2*g*s*sp.sin(al)); lam1=h/(m*v1)
dphi_exact=2*sp.pi*l*(1/lam - 1/lam1)   # path ABC: AB at lam0, path ADC: DC at lam1 (AD,BC identical)
dphi_exact=dphi_exact.subs(v0,h/(m*lam))
ser=sp.series(dphi_exact.subs(g,sp.Symbol('eps')*g),sp.Symbol('eps'),0,2).removeO().subs(sp.Symbol('eps'),1)
print('first order dphi =',sp.simplify(ser))

# ---- new (4): general contamination fraction w ----
ph,w=sp.symbols('phi w',positive=True)
Iw=(1-w)*(1+sp.cos(ph))+w*(1+sp.cos(ph/2))
dI=sp.diff(Iw,ph)
fac=-sp.sin(ph/2)*(2*(1-w)*sp.cos(ph/2)+w/2)
print('dI/dphi - factored form =',sp.simplify(sp.expand_trig(dI-fac)))
cst=-w/(4*(1-w))
Imin=sp.simplify((2*(1-w)*sp.Symbol('c')**2+w*sp.Symbol('c')+w).subs(sp.Symbol('c'),cst))
print('I at min =',sp.simplify(Imin),'  ratio formula check:',sp.simplify(Imin/2-(w-w**2/(8*(1-w)))/2))
print('w=1/2 ratio:',sp.nsimplify(Imin.subs(w,sp.Rational(1,2))/2),'  boundary |c|=1 at w=',sp.solve(sp.Eq(w/(4*(1-w)),1),w))
for wv in [0.05,0.3,0.5,0.75,0.85]:
    ps=np.linspace(1e-6,2*np.pi,400001); Iv=(1-wv)*(1+np.cos(ps))+wv*(1+np.cos(ps/2))
    dd=np.diff(Iv); k=np.argmax(dd>0) if (dd>0).any() else len(ps)-1
    pred=2*np.arccos(-wv/(4*(1-wv))) if wv<0.8 else 2*np.pi
    print(' w=%.2f first min phi num %.5f pred %.5f ratio num %.5f   lambda\'/lambda0 = %.4f'%(wv,ps[k],pred,Iv[k]/2,np.pi/pred))
