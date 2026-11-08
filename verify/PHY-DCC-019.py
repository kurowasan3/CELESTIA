# S0511-1 verification: general ladder (shunt a, series b, termination t), 3 nodes, nodal analysis
import sympy as sp
a,b,t,E,RL=sp.symbols('a b t E R_L',positive=True)
def vout(U,term=None,bb=None,load=None):
    term = t if term is None else term; bb = b if bb is None else bb
    V=sp.symbols('V0:3')
    eqs=[V[0]/term+(V[0]-U[0])/a+(V[0]-V[1])/bb,
         (V[1]-U[1])/a+(V[1]-V[0])/bb+(V[1]-V[2])/bb,
         (V[2]-U[2])/a+(V[2]-V[1])/bb+(V[2]/load if load else 0)]
    s=sp.solve(eqs,V,dict=True)[0]; return sp.simplify(s[V[2]])
# (1) self-similarity: b + a*t/(a+t) = t
bsol=sp.solve(sp.Eq(b+a*t/(a+t),t),b)[0]; print('b =',sp.factor(bsol))
U=sp.symbols('U')
w=[sp.factor(vout([U if i==j else 0 for i in range(3)],bb=bsol)) for j in range(3)]
print('weights P0,P1,P2:',w)
assert sp.simplify(w[0]-U*t*a**2/(a+t)**3)==0 and sp.simplify(w[1]-U*t*a/(a+t)**2)==0 and sp.simplify(w[2]-U*t/(a+t))==0
# (3) ternary: a/(a+t)=1/3 -> t=2a
tt=sp.solve(sp.Eq(a/(a+t),sp.Rational(1,3)),t)[0]; bt=sp.factor(bsol.subs(t,tt)); print('t=',tt,' b=',bt)
d=sp.symbols('d0:3')
V=vout([d[0]*E/2,d[1]*E/2,d[2]*E/2],term=tt,bb=bt)
print('Vout =',sp.factor(V)); assert sp.simplify(V-E*(9*d[2]+3*d[1]+d[0])/27)==0
# (4) output resistance: all switches 0, inject test current -> V/I
I=sp.symbols('I'); Vs=sp.symbols('W0:3')
eqs=[Vs[0]/tt+Vs[0]/a+(Vs[0]-Vs[1])/bt,Vs[1]/a+(Vs[1]-Vs[0])/bt+(Vs[1]-Vs[2])/bt,Vs[2]/a+(Vs[2]-Vs[1])/bt-I]
s=sp.solve(eqs,Vs,dict=True)[0]; Rout=sp.simplify(s[Vs[2]]/I); print('Rout=',Rout)
assert sp.simplify(Rout-2*a/3)==0
# load check by direct nodal analysis
VL=vout([E,E,E],term=tt,bb=bt,load=RL); print('VL(222)=',sp.factor(VL))
assert sp.simplify(VL-sp.Rational(26,27)*E*RL/(RL+2*a/3))==0
# (5) error < half LSB (E/54): 26/27 E *(2a/3)/(RL+2a/3) < E/54  -> RL > 34 a
x=sp.symbols('x',positive=True)
print('RL threshold:',sp.solve(sp.Eq(sp.Rational(26,27)*(2*a/3)/(x+2*a/3),sp.Rational(1,54)),x))
# binary check: t=a, b=a/2 reproduces R-2R
print('binary b:',bsol.subs(t,a))
print('ALL OK')
