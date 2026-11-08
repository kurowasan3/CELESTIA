# S0315-1: bore pushed by a paddle in a channel (shallow-water wave speed by momentum method)
import sympy as sp, numpy as np
rho,g,b,h,d,v,U,t=sp.symbols('rho g b h Delta v U t',positive=True)
H=h+d
mass=sp.Eq(H*(U-v),h*U)                      # (1)
vsol=sp.solve(mass,v)[0]; print('(1) v =',vsol)
# (2) momentum: d/dt[rho b h U t v] = 1/2 rho g b (H^2-h^2)
mom=sp.Eq(rho*b*h*U*v, sp.Rational(1,2)*rho*g*b*(H**2-h**2))
U2=sp.solve(mom.subs(v,vsol),U)
print('(2) U =',[sp.factor(u) for u in U2])
Uex=sp.sqrt(g*H*(2*h+d)/(2*h))
print('    check:',sp.simplify(U2[0]-Uex))
print('    small-amp series U:',sp.series(Uex,d,0,2))
# (3) energy budget
F=sp.Rational(1,2)*rho*g*b*H**2   # paddle force (the far fixed boundary does no work)
W=F*v
KE=sp.Rational(1,2)*rho*b*h*U*v**2
PE=sp.Rational(1,2)*rho*g*b*(H**2*(U-v)-h**2*U)
D=sp.simplify((W-KE-PE).subs(v,vsol))
print('(3) dissipation rate =',sp.factor(D))
Dsub=sp.factor(sp.simplify(D.subs(U,Uex)))
print('    with U:',Dsub)
target=rho*g*b*U*d**3/(4*H)
print('    target rho g b U Delta^3/(4(h+Delta)) diff:',sp.simplify(D-target))
ratio=sp.simplify((D/W).subs(v,vsol))
print('    D/W =',sp.factor(ratio), ' series:',sp.series(ratio,d,0,3))
# fraction of work to KE and PE
print('    KE/W=',sp.factor(sp.simplify((KE/W).subs(v,vsol))),' PE/W=',sp.factor(sp.simplify((PE/W).subs(v,vsol))))
# numeric: h=1 m, Delta=0.2 m
num={g:9.8,h:1.0,d:0.2}
print('    numeric U=%.4f sqrt(gh)=%.4f v=%.4f D/W=%.4f'%(float(Uex.subs(num)),np.sqrt(9.8),float(vsol.subs(U,Uex).subs(num)),float(ratio.subs(v,vsol).subs(U,Uex).subs(num))))
# (4) energy-conservation (wrong) speed
UE=sp.solve(sp.Eq(D,0),U); print('(4) U_E candidates:',UE)
UEx=[u for u in UE if sp.simplify(u-H*sp.sqrt(g/h))==0]
assert UEx, UE
UEs=H*sp.sqrt(g/h)
# (4a) U_E^2/U^2 = 2(h+Delta)/(2h+Delta) = 1 + Delta/(2h+Delta) > 1
r2=sp.simplify(UEs**2/Uex**2); print('(4a) U_E^2/U^2 =',sp.factor(r2))
assert sp.simplify(r2-(1+d/(2*h+d)))==0
# (4b) U_E/U = 1 + Delta/(4h) + O(Delta^2)
eps_=sp.symbols('epsilon',positive=True)
ser=sp.series(sp.sqrt(r2).subs(d,eps_*h),eps_,0,2).removeO(); print('(4b) U_E/U =',ser)
assert sp.simplify(ser-(1+eps_/4))==0
print('    U/sqrt(gh) =',sp.series(Uex.subs(d,eps_*h)/sp.sqrt(g*h),eps_,0,2),'  U_E/sqrt(gh) =',sp.expand(UEs.subs(d,eps_*h)/sp.sqrt(g*h)))
# (4c) 各項の次数：P,G は Delta の1次，P-G と K は2次，Q/(P-G)=Delta/(2(h+Delta))
Pv=sp.simplify(W.subs(v,vsol)); Gv=sp.simplify(PE.subs(v,vsol))
assert sp.simplify(Pv-sp.Rational(1,2)*rho*g*b*H*d*U)==0 and sp.simplify(Gv-sp.Rational(1,2)*rho*g*b*h*d*U)==0
assert sp.simplify((Pv-Gv)-sp.Rational(1,2)*rho*g*b*d**2*U)==0
Kv=sp.simplify(KE.subs(v,vsol)); assert sp.simplify(Kv-sp.Rational(1,2)*rho*b*h*d**2*U**3/H**2)==0
Dv=sp.simplify(D.subs(U,Uex))
qr=sp.simplify(Dv/(Pv-Gv).subs(U,Uex)); print('    Q/(P-G) =',qr)
assert sp.simplify(qr-d/(2*H))==0
qp=sp.simplify(Dv/Pv.subs(U,Uex)); print('    Q/P =',qp); assert sp.simplify(qp-d**2/(2*H**2))==0
# 数値：U_E/U と 1+Delta/(4h)
for dv in [0.01,0.05,0.2]:
    ex=float((UEs/Uex).subs({g:9.8,h:1.0,d:dv})); print('    Delta/h=',dv,' U_E/U=',ex,' 1+Delta/4h=',1+dv/4)
    assert abs(ex-(1+dv/4))<dv**2
    assert ex>1
print('PHY-WAV-020 (4) OK')
