# N62-1 recoil shift, resonance absorption, moving source cone
import sympy as sp, numpy as np
M,c,dE,E,u,th=sp.symbols('M c DeltaE E u theta',positive=True)
# (1) emission at rest: E + E^2/(2Mc^2) = dE
Ee=sp.solve(sp.Eq(E+E**2/(2*M*c**2),dE),E); Ee=[s for s in Ee if s.subs({M:1,c:1,dE:sp.Rational(1,10)})>0][0]
eps=sp.Symbol('epsilon',positive=True)
print('(1) Ee series:',sp.series(Ee.subs(dE,eps*M*c**2)/(M*c**2),eps,0,3))
# (2) absorption at rest: dE + E^2/(2Mc^2) = E
Ea=sp.solve(sp.Eq(dE+E**2/(2*M*c**2),E),E); Ea=min(Ea,key=lambda s: s.subs({M:1,c:1,dE:sp.Rational(1,10)}))
print('(2) Ea series:',sp.series(Ea.subs(dE,eps*M*c**2)/(M*c**2),eps,0,3))
print('    Ea-Ee series:',sp.series((Ea-Ee).subs(dE,eps*M*c**2)/(M*c**2),eps,0,3))
# (3) moving source: momentum/energy -> dE = E(1 - u cos/c + E/(2Mc^2))
ux,uy=sp.symbols('ux uy')
up2=(u-E*sp.cos(th)/(M*c))**2+(E*sp.sin(th)/(M*c))**2
rel=sp.simplify(M*u**2/2+dE-(M*up2/2+E))
print('(3) energy balance -> dE =',sp.solve(rel,dE)[0].expand())
cos_res=sp.solve(sp.Eq(sp.solve(rel,dE)[0], sp.solve(sp.Eq(dE+E**2/(2*M*c**2),E),dE)[0]),sp.cos(th))
print('    cos(theta) for resonance:',cos_res)
# (4) numbers 57Fe
amu=1.66e-27; Mv=57*amu; cv=3.0e8; ev=1.6e-19; dEv=14.4e3*ev
Mc2=Mv*cv**2
Eev=float(Ee.subs({M:Mv,c:cv,dE:dEv})); Eav=float(Ea.subs({M:Mv,c:cv,dE:dEv}))
print('(4) Mc2 [eV]',Mc2/ev,' Ea-Ee [eV]',(Eav-Eev)/ev,' approx',dEv**2/Mc2/ev)
umin=Eav/(Mv*cv); print('    umin',umin,' theta at u=100:',np.degrees(np.arccos(umin/100)))
# direct check by numerically solving conservation for u=100, theta found -> photon energy equals Ea
thv=np.arccos(umin/100); from scipy.optimize import brentq
f=lambda Ev: Mv*100**2/2+dEv-(Mv*((100-Ev*np.cos(thv)/(Mv*cv))**2+(Ev*np.sin(thv)/(Mv*cv))**2)/2+Ev)
Esol=ev*brentq(lambda x: f(x*ev)/ev,0.9*dEv/ev,1.1*dEv/ev,xtol=1e-12); print('    emitted E - Ea [eV]:',(Esol-Eav)/ev)
# (5) hydrogen
MH=1.7e-27; dEH=1.9*ev; print('(5) umin H',dEH/(MH*cv), ' gap H [eV]',dEH**2/(MH*cv**2)/ev)
print('    vrms H 300K', np.sqrt(3*8.3*300/1.0e-3))
