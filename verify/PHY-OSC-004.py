# M29-1: block M on smooth floor, top of its smooth surface is a circular arc radius R (center fixed to block), top at height H.
# ball m launched with v0 along floor toward block (block at rest). At top: velocities, N, conditions.
import sympy as sp
m,M,g,R,H,v0=sp.symbols('m M g R H v_0',positive=True)
vb,V=sp.symbols('v_b V')
sols=sp.solve([sp.Eq(m*v0,m*vb+M*V), sp.Eq(m*v0**2/2, m*vb**2/2+M*V**2/2+m*g*H)],[vb,V],dict=True)
for s in sols: print({k:sp.simplify(v) for k,v in s.items()}, ' u^2=',sp.simplify((s[vb]-s[V])**2))
u2=v0**2-2*g*H*(M+m)/M
assert all(sp.simplify((s[vb]-s[V])**2-u2)==0 for s in sols)
# Lagrangian on the arc: block position X, ball at angle phi from top (on arc, center C fixed in block at height H-R)
t=sp.symbols('t'); X=sp.Function('X')(t); ph=sp.Function('phi')(t)
xb=X+R*sp.sin(ph); yb=H-R+R*sp.cos(ph)
L=sp.Rational(1,2)*M*X.diff(t)**2+sp.Rational(1,2)*m*(xb.diff(t)**2+yb.diff(t)**2)-m*g*yb
eqs=[sp.diff(L.diff(v.diff(t)),t)-L.diff(v) for v in (X,ph)]
Xdd,phdd=sp.symbols('Xdd phdd'); Xd,phd=sp.symbols('Xd phd'); p=sp.symbols('p')
rep={X.diff(t,2):Xdd, ph.diff(t,2):phdd}
eqs=[e.subs(rep) for e in eqs]
rep2={X.diff(t):Xd, ph.diff(t):phd, ph:p}
eqs=[sp.simplify(e.subs(rep2)) for e in eqs]
sol=sp.solve(eqs,[Xdd,phdd],dict=True)[0]
# normal force: from block equation M Xdd = N sin(phi) (ball pushes block along -normal... ) -> N = M Xdd / sin(phi); use ball radial eq instead:
# ball acceleration a = (Xdd,0) + relative accel; N along outward normal n=(sin p, cos p): m a = N n + (0,-mg)
ax=Xdd+sp.diff(R*sp.sin(ph),t,2).subs(rep).subs(rep2)
ay=sp.diff(R*sp.cos(ph),t,2).subs(rep).subs(rep2)
Nexpr=sp.simplify((m*ax*sp.sin(p)+(m*ay+m*g)*sp.cos(p)).subs(sol))
N_top=sp.simplify(Nexpr.subs(p,0)); print('N at top =',N_top, '  Xdd at top=',sp.simplify(sol[Xdd].subs(p,0)))
uu=sp.symbols('u')
assert sp.simplify(N_top.subs(phd,uu/R)-m*(g-uu**2/R))==0
assert sp.simplify(sol[Xdd].subs(p,0))==0
# conditions
print('reach top: v0^2 >',sp.simplify(2*g*H*(M+m)/M),'; not leave: v0^2 <=',sp.simplify(g*R+2*g*H*(M+m)/M))
# numeric case M=m, H=R
d={M:m,H:R}
print('M=m,H=R: v0^2 in (',sp.simplify((2*g*H*(M+m)/M).subs(d)),',',sp.simplify((g*R+2*g*H*(M+m)/M).subs(d)),']')
# example: at v0^2=4.5gR, N=?
print('v0^2=9gR/2 -> N=',sp.simplify((m*(g-u2/R)).subs(d).subs(v0,sp.sqrt(sp.Rational(9,2)*g*R))))
print('OK')
