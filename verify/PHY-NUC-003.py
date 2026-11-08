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
# ---- new (5): after-process velocities are perpendicular, KE conserved (equal-mass elastic collision) ----
Ea_s,cth=sp.symbols('E_a costh',positive=True)
umin_s=Ea_s/(M*c)
cos_cond=Ea_s/(M*u*c)
n=sp.Matrix([cos_cond, sp.sqrt(1-cos_cond**2)])
U=sp.Matrix([u,0])
vB=umin_s*n
vA=U-umin_s*n
print('vA.vB =',sp.simplify(vA.dot(vB)))
print('|vA|^2 - (u^2-umin^2) =',sp.simplify(vA.dot(vA)-(u**2-umin_s**2)))
print('KE balance =',sp.simplify(M*u**2/2 - M*vA.dot(vA)/2 - M*vB.dot(vB)/2))
# full energy check incl. internal energy: dE(Ea) from absorption relation
dE_of_Ea=Ea_s-Ea_s**2/(2*M*c**2)
print('emission energy balance (A):',sp.simplify(M*u**2/2+dE_of_Ea-(M*vA.dot(vA)/2+Ea_s)))
print('absorption energy balance (B):',sp.simplify(Ea_s-(dE_of_Ea+M*vB.dot(vB)/2)))
print('u=umin -> vA =',sp.simplify(vA.subs(u,umin_s)).T)
# numeric spot check with arbitrary values
vals={M:57*1.66e-27,c:3e8,Ea_s:14.4e3*1.6e-19,u:100.0}
print('numeric vA.vB',float(vA.dot(vB).subs(vals)),' |vA|',float(sp.sqrt(vA.dot(vA)).subs(vals)))
