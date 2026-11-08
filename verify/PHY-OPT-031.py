# J094-1: Michelson interferometer, tilted glass plate: OPL excess t(sqrt(n^2-sin^2 i)-cos i)
import sympy as sp, numpy as np
n,t,i,lam=sp.symbols('n t i lambda',positive=True)
r=sp.asin(sp.sin(i)/n)
# geometric derivation: ray in plate length t/cos r, compare with straight air path between perpendicular planes
# entry point (0,0), plate faces normal to x' axis rotated by i. Use coordinates: beam along X. Plate normal at angle i to X.
# exit point in plate: travel t/cos r along direction at angle (i - r)?? -> beam direction inside makes angle (i-r) with X
Lin=t/sp.cos(r)
proj=Lin*sp.cos(i-r)   # X-distance covered inside plate
opl_with=n*Lin           # then air path continues; compare at a plane X=const beyond: without plate path = proj
excess=sp.simplify(opl_with-proj)
print('excess - formula:', sp.simplify(sp.expand_trig(excess-(t*(sp.sqrt(n**2-sp.sin(i)**2)-sp.cos(i))))))
# numeric check of above identity
f=lambda nv,iv: nv/np.cos(np.arcsin(np.sin(iv)/nv))-np.cos(iv-np.arcsin(np.sin(iv)/nv))/np.cos(np.arcsin(np.sin(iv)/nv))
for nv,iv in [(1.5,0.3),(1.7,1.0)]:
    print(f(nv,iv), np.sqrt(nv**2-np.sin(iv)**2)-np.cos(iv))
# Independent check via wavefront phase: plane wave k along X, phase at a far plane X=X0 for path through plate.
# Exit ray parallel to X, so OPL to plane X0 = n*Lin + (X0 - proj); without plate = X0 -> excess same. Also check lateral shift irrelevant (plane wave).
N=2*t*((sp.sqrt(n**2-sp.sin(i)**2)-sp.cos(i))-(n-1))/lam
print('small angle N ~', sp.simplify(sp.series(N,i,0,4).removeO()))
# numeric: t=1.00 mm, lam=500 nm, i=30deg
tv,lv=1.0e-3,5.0e-7
for nv in [1.5,1.52]:
    print('n',nv,'N(30deg)=',float(N.subs({n:nv,t:tv,lam:lv,i:sp.pi/6})))
# inverse: given N=200 at 30 deg
Nv=200; K=Nv*lv/(2*tv); s=0.5; c=np.cos(np.pi/6)
cp=c-1+K
nsol=(s**2+cp**2)/(-2*cp)
print('K',K,'c-1+K',cp,'n =',nsol)
print('check N(nsol)=',float(N.subs({n:nsol,t:tv,lam:lv,i:sp.pi/6})))
# symbolic inverse
Ks=sp.symbols('K'); ns=sp.symbols('n_s')
sol=sp.solve(sp.Eq(sp.sqrt(ns**2-sp.sin(i)**2), ns-1+sp.cos(i)+Ks),ns)
print('inverse sym:',[sp.simplify(s_) for s_ in sol])
# small-angle approximate n from N at 30 deg (to show inaccuracy)
ia=np.pi/6; napprox=tv*ia**2/(tv*ia**2-Nv*lv); print('small-angle approx n =',napprox)
# mirror compensation: OPL change of round trip = N lam -> mirror moves N lam/2 toward H
print('mirror shift (mm)', Nv*lv/2*1e3)
