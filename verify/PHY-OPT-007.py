# N06-1: oil film (1.40) on water (1.33), reflected white light
import sympy as sp, numpy as np
n=sp.Rational(14,10)
m=sp.symbols('m')
ms=sp.solve(sp.Eq((m+sp.Rational(1,2))*640,(m+sp.Rational(3,2))*457),m); print("m approx",[float(v) for v in ms])
d=sp.Rational(5,2)*640/(2*n); print("d nm =",float(d))
print("maxima nm:",[float(2*n*d/(k+sp.Rational(1,2))) for k in range(0,6)])
print("minima nm:",[float(2*n*d/k) for k in range(1,6)])
# thickness ranges with no maximum in [400,700]
ds=np.linspace(1,400,399001); ok=[]
for dv in ds:
    lams=[2*1.4*dv/(k+0.5) for k in range(0,20)]
    if not any(400<=L<=700 for L in lams): ok.append(dv)
ok=np.array(ok); br=np.where(np.diff(ok)>0.002)[0]
print("no-max ranges:", ok[0], ok[br], ok[br+1], ok[-1])
print("bounds:",400/5.6,700/5.6,1200/5.6)
# oblique 60 deg
lam60=2*float(d)*np.sqrt(1.96-np.sin(np.radians(60))**2)/2.5; print("lambda at 60deg:",lam60)
# exact two-beam reflectance check with Fresnel amplitudes (normal incidence) -> max positions
r1=(1-1.4)/(1+1.4); r2=(1.4-1.33)/(1.4+1.33)
L=np.linspace(380,1200,82001); ph=4*np.pi*1.4*float(d)/L
R=np.abs(r1+r2*np.exp(1j*ph))**2
pk=L[1:-1][(R[1:-1]>R[:-2])&(R[1:-1]>R[2:])]; print("Fresnel maxima:",pk)
