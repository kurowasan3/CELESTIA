# Q173-1 rotating thick metal cylinder in axial B
import sympy as sp
r,a,b,e,m,B,eps,w,h=sp.symbols('r a b e m B epsilon0 omega h',real=True)
# electron (charge -e) co-rotating, omega vector parallel to B (z). radial eq: -m w^2 r = (-e)E_r + (-e)(v x B)_r, (v x B)_r = w r B
Er=sp.symbols('E_r')
E=sp.solve(sp.Eq(-m*w**2*r, -e*Er - e*w*r*B),Er)[0]; print('E_r=',sp.expand(E))
rho=sp.simplify(eps*sp.diff(r*E,r)/r); print('rho=',sp.factor(rho))
sa=sp.simplify(eps*E.subs(r,a)); sb=sp.simplify(-eps*E.subs(r,b)); print('sigma_a=',sp.factor(sa),' sigma_b=',sp.factor(sb))
print('neutral:',sp.simplify(sa*2*sp.pi*a+rho*sp.pi*(b**2-a**2)+sb*2*sp.pi*b))
# Gauss check: flux through cylinder radius r = enclosed / eps
enc=sa*2*sp.pi*a*h+rho*sp.pi*(r**2-a**2)*h
print('gauss:',sp.simplify(E*2*sp.pi*r*h-enc/eps))
dV=sp.simplify(sp.integrate(E,(r,a,b))); print('V(a)-V(b)=',sp.factor(dV))
print('zero rho at w=',sp.solve(sp.Eq(rho,0),w))
# small-omega limit (m->0)
print('m->0: rho',sp.factor(rho.subs(m,0)),' sa',sp.factor(sa.subs(m,0)),' sb',sp.factor(sb.subs(m,0)),' dV',sp.factor(dV.subs(m,0)))
# numbers: B=1T
print('wc=',1.602e-19*1.0/9.11e-31)
# numeric example rho magnitude for w=100 rad/s, B=1 T
print('rho num', 2*8.85e-12*100*1.0, 'electrons per m^3',2*8.85e-12*100/1.602e-19)
