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
print('    U_E numeric:',[float(u.subs(num)) for u in UE if u.subs(num).is_positive])
