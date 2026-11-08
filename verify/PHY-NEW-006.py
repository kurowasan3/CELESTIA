# M22-1: block P (m) on smooth slope of wedge Q (M) on a rough floor (mu0 static, mu kinetic)
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,M,g,th,mu,A,N=sp.symbols('m M g theta mu A N',positive=True)
# (1)(2) wedge at rest
N0=m*g*sp.cos(th); R0=M*g+N0*sp.cos(th); f0=N0*sp.sin(th)
print('(1) R0=',sp.simplify(R0),' f0=',sp.simplify(f0),' (2) mu0 >=',sp.simplify(f0/R0))
# (3) inertial-force method in Q frame (Q accelerates to the left with A)
sol=sp.solve([sp.Eq(N+m*A*sp.sin(th),m*g*sp.cos(th)), sp.Eq(M*A,N*sp.sin(th)-mu*(M*g+N*sp.cos(th)))],[A,N],dict=True)[0]
vals={m:1,M:1,th:sp.pi/4,mu:sp.Rational(1,5)}
print('(3) inertial method A=',sp.nsimplify(sp.simplify(sol[A].subs(vals))),' N=',sp.nsimplify(sp.simplify(sol[N].subs(vals))))
print('    general A=',sp.simplify(sol[A]))
# independent: ground frame, unknowns ax, ay (P), Xdd (Q, +x right), N, R
ax,ay,Xdd,R=sp.symbols('a_x a_y Xdd R')
n=sp.Matrix([sp.sin(th),sp.cos(th)])       # outward normal of slope (slope descends to the right)
d=sp.Matrix([sp.cos(th),-sp.sin(th)])      # down-slope direction
eqs=[sp.Eq(m*ax,N*n[0]), sp.Eq(m*ay,N*n[1]-m*g),
     sp.Eq(M*Xdd,-N*n[0]+mu*R),            # Q moves left (Xdd<0) so kinetic friction to the right
     sp.Eq(R,M*g+N*n[1]),
     sp.Eq((sp.Matrix([ax-Xdd,ay]).T*n)[0],0)]  # relative acceleration tangent to slope
s2=sp.solve(eqs,[ax,ay,Xdd,N,R],dict=True)[0]
print('ground-frame: Xdd=',sp.nsimplify(sp.simplify(s2[Xdd].subs(vals))),' N=',sp.nsimplify(sp.simplify(s2[N].subs(vals))),' R=',sp.nsimplify(sp.simplify(s2[R].subs(vals))))
arel=sp.simplify(((sp.Matrix([s2[ax]-s2[Xdd],s2[ay]]).T*d)[0]).subs(vals))
print('(4) a_rel =',sp.nsimplify(arel),' expected 4*sqrt(2)*g/7')
L=sp.symbols('L',positive=True)
print('(4) time to slide L:',sp.simplify(sp.sqrt(2*L/arel)))
print('static threshold numeric', sp.simplify((f0/R0).subs(vals)))
# kinetic threshold for A>0
print('A>0 iff mu <', sp.simplify(sp.solve(sp.Eq(sol[A],0),mu)[0]))
