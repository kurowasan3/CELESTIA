# M23-1 verification: block A (mass M) with pulley at right top corner, B (m) on top, C (m) in vertical groove on A's right face.
# floor friction coefficient mu (static = kinetic). A moves left with acceleration al (>0 to the left).
import sympy as sp
M,m,g,mu,h=sp.symbols('M m g mu h',positive=True)
al,b,T,NC,R=sp.symbols('alpha b T N_C R')
# Method 1: ground frame, x to the right positive.
# A accel: -al ; B accel x: -al + b ; C accel: x: -al, y: -b (down)
eqs=[
 sp.Eq(m*(-al+b), T),                 # B horizontal: tension toward pulley (right)
 sp.Eq(m*(-b), T-m*g),                # C vertical
 sp.Eq(m*(-al), -NC),                 # C horizontal: groove wall pushes C left with NC
 sp.Eq(M*(-al), -T + NC + mu*R),      # A horizontal: pulley pulled left by T, C pushes A right with NC, friction right (A moves left)
 sp.Eq(R, (M+m)*g + T),               # A vertical: weight, B normal mg, pulley vertical string T
]
s=sp.solve(eqs,[al,b,T,NC,R],dict=True)[0]
alpha=sp.factor(s[al]); print('alpha=',alpha)
print('T=',sp.factor(s[T]),' b=',sp.factor(s[b]),' R=',sp.factor(s[R]))
ans_alpha = g*(m-mu*(2*M+3*m))/(2*M+3*m-mu*m)
assert sp.simplify(alpha-ans_alpha)==0
# threshold: A static
s0=sp.solve([e.subs(al,0) for e in eqs[:3]]+[sp.Eq(R,(M+m)*g+T)],[b,T,NC,R],dict=True)[0]
fneed=s0[T]  # friction needed = T (NC=0)
print('static: T=',s0[T],' R=',s0[R],' mu_min=',sp.simplify(fneed/s0[R]))
assert sp.simplify(fneed/s0[R]-m/(2*M+3*m))==0
# Method 2: system total vertical: R_total = (M+2m)g - m b
assert sp.simplify(s[R]-((M+2*m)*g-m*s[b]))==0
# (3) frictionless: displacement of A while C descends h  = alpha*h/b
dA=sp.simplify((s[al]/s[b]*h).subs(mu,0)); print('dA (mu=0)=',dA)
assert sp.simplify(dA-m*h/(M+2*m))==0
# (2) answers for given mu
print('T(mu)=',sp.factor(s[T]))
ansT=m*g*(1+mu)*(M+m)/(2*M+3*m-mu*m)
assert sp.simplify(s[T]-ansT)==0
# physical checks: for 0<=mu<m/(2M+3m): alpha>0, T>0, b>0
import random
for _ in range(2000):
    Mv=random.uniform(0.1,10); mv=random.uniform(0.1,10); muv=random.uniform(0,mv/(2*Mv+3*mv))
    d={M:Mv,m:mv,mu:muv,g:9.8}
    assert s[al].subs(d)>=-1e-12 and s[T].subs(d)>0 and s[b].subs(d)>0 and s[R].subs(d)>0
# (4) numeric: M=2m, mu=1/14 ->
d={M:2*m,mu:sp.Rational(1,14)}
print('M=2m,mu=1/14 alpha=',sp.simplify(s[al].subs(d)),'b=',sp.simplify(s[b].subs(d)),'dA=',sp.simplify((s[al]/s[b]*h).subs(d)))
print('OK')
