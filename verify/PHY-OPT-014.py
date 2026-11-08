# P49-1: two incoherent source slits (separation s) -> double slit (d) -> screen; exact geometry numeric
import numpy as np, sympy as sp
lam,L,R=600e-9,2.0,1.0
def I_single(ys,d,x,extra1=0.0):
    # ys: source y; slits at +-d/2; path via S1 + extra1 (glass in front of S1)
    p1=np.hypot(L,ys-d/2)+np.hypot(R,x-d/2)+extra1
    p2=np.hypot(L,ys+d/2)+np.hypot(R,x+d/2)
    return np.abs(np.exp(2j*np.pi*p1/lam)+np.exp(2j*np.pi*p2/lam))**2
def vis(I): return (I.max()-I.min())/(I.max()+I.min())
s=1.0e-3
x=np.linspace(-5e-3,5e-3,20001)
# (1) single source at +s/2: zeroth-order fringe near x=-sR/(2L)
d=0.3e-3
I=I_single(s/2,d,x); xs=x[np.abs(x)<1e-3]; Is=I[np.abs(x)<1e-3]
print('zeroth fringe near', xs[np.argmax(Is)], 'expected', -s*R/(2*L), 'fringe spacing', lam*R/d)
# (2)(3) visibility vs d for two sources
ds=np.linspace(0.05e-3,1.5e-3,300)
V=[vis(I_single(s/2,dd,x)+I_single(-s/2,dd,x)) for dd in ds]
i0=np.argmax(np.array(V)<0.02); print('first vanish d ~',ds[i0],'expected',lam*L/(2*s))
print('vis at d=lam L/s', vis(I_single(s/2,lam*L/s,x)+I_single(-s/2,lam*L/s,x)))
# inverse with d1=0.60mm
print('s from d1=0.6mm:', lam*L/(2*0.60e-3))
# (4) glass plate in front of S1: extra path (n-1)t
dd=0.3e-3; ex=(1.5-1)*2e-6
Ia=I_single(s/2,dd,x)+I_single(-s/2,dd,x); Ib=I_single(s/2,dd,x,ex)+I_single(-s/2,dd,x,ex)
print('vis without/with plate at S1', vis(Ia), vis(Ib))
k=np.argmax(Ia*(np.abs(x)<lam*R/dd/2)); k2=np.argmax(Ib*(np.abs(x-ex*R/dd)<lam*R/dd/2))
print('shift of central max', x[k2]-x[k], 'expected +(n-1)tR/d =', ex*R/dd)
# plate in front of source A1: adds same phase to both paths -> identical pattern (trivially), check
# (5) continuous source width w
w=1.0e-3
def Icont(d):
    ys=np.linspace(-w/2,w/2,401); return sum(I_single(y,d,x) for y in ys)
Vc=[(dd,vis(Icont(dd))) for dd in np.linspace(0.6e-3,1.6e-3,21)]
print('continuous vis min at', min(Vc,key=lambda t:t[1]), 'expected d=lam L/w', lam*L/w)
# symbolic: path difference
y,ds_,Ls,x_,Rs=sp.symbols('y d L x R',positive=True)
diff=sp.series(sp.sqrt(Ls**2+(y+ds_/2)**2)-sp.sqrt(Ls**2+(y-ds_/2)**2),Ls,sp.oo,2).removeO()
print('A->S2 minus A->S1 ~',sp.simplify(diff))
