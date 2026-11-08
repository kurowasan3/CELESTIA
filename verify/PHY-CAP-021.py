# Two isolated parallel plates with charges Q1 (lower), Q2 (upper); partial dielectric slab thickness t on lower plate
import sympy as sp
Q1,Q2,e0,er,L,M,d,t,b,m=sp.symbols('Q1 Q2 epsilon0 epsilon_r L M d t b m',positive=True)
S=L*M
# (1) sheet superposition: unknown surface charges a(outer low), c(inner low), e(inner up), f(outer up)
a,c,e,f=sp.symbols('a c e f')
# field inside lower conductor (z<0..0) and upper conductor =0 : field = sum sign*sigma/(2e0)
eqs=[sp.Eq(a+c,Q1),sp.Eq(e+f,Q2),
     sp.Eq((-a + c + e + f),0),   # inside lower conductor: sheets above push down: field(+z) = (a - c - e - f)/(2 e0 S)=0
     sp.Eq((a + c + e - f),0)]    # inside upper conductor
print('(1)',sp.solve(eqs,[a,c,e,f]))
Q=(Q1-Q2)/2
dp=d-t+t/er
sa,sb=sp.symbols('sa sb')
sol=sp.solve([sp.Eq(sa*dp,sb*d),sp.Eq(sa*M*b+sb*M*(L-b),Q)],[sa,sb])
print('(2) sa=',sp.simplify(sol[sa]),' sb=',sp.simplify(sol[sb]))
C=e0*M*(b/dp+(L-b)/d)
V=sp.simplify(Q/C); print('(3) V=',V, ' check V=sa*dp/e0:',sp.simplify(V-sol[sa]*dp/e0))
print('    bound charge top surface density:',sp.simplify(sol[sa]*(1-1/er)))
U=Q**2/(2*C)
F=sp.simplify(-sp.diff(U,b)); print('(4) F=',sp.factor(F))
Fexp=Q**2*e0*M*(1/dp-1/d)/(2*C**2); print('    matches formula:',sp.simplify(F-Fexp)==0)
# (5) speed when fully inserted from b=0 (frictionless), slab mass m
v=sp.sqrt(2*(U.subs(b,0)-U.subs(b,L))/m); print('(5) v^2=',sp.factor(sp.simplify(v**2)))
# numeric sanity: Q1=Q2 -> no force
print('Q1=Q2 force:',sp.simplify(F.subs(Q2,Q1)))
# independent check: force from field energy density integration (vacuum region energy + dielectric region energy)
Ea=sol[sa]/e0; Eb=sol[sb]/e0
Ufield=sp.Rational(1,2)*e0*Eb**2*M*(L-b)*d + M*b*(sp.Rational(1,2)*e0*Ea**2*(d-t)+sp.Rational(1,2)*e0*er*(Ea/er)**2*t)
print('field-energy route equal U:',sp.simplify(Ufield-U)==0)
