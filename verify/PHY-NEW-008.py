# M24-1: A(m), B(2m) on string1 over movable massless pulley P; P and C(M) on string2 over fixed pulley Q.
import sympy as sp
m,M,g=sp.symbols('m M g',positive=True)
aA,aB,aP,T1,T2,be,al=sp.symbols('a_A a_B a_P T1 T2 beta alpha')
# Method 1 (ground frame, up positive, constraint)
s=sp.solve([m*aA-(T1-m*g),2*m*aB-(T1-2*m*g),aA+aB-2*aP,M*(-aP)-(T2-M*g),T2-2*T1],[aA,aB,aP,T1,T2],dict=True)[0]
for k in s: print(k,sp.factor(s[k]))
# Method 2 (P frame, inertial force -m*aP; A rel accel up beta, B rel down beta)
s2=sp.solve([m*be-(T1-m*g-m*al), 2*m*be-(2*m*g+2*m*al-T1), M*al-(M*g-2*T1)],[be,al,T1],dict=True)[0]
# al = C's downward accel = P's upward accel
assert sp.simplify(s2[al]-s[aP])==0 and sp.simplify(s2[T1]-s[T1])==0
assert sp.simplify(s[aA]-(s2[al]+s2[be]))==0 and sp.simplify(s[aB]-(s2[al]-s2[be]))==0
# thresholds
for k in [aA,aP,aB]:
    print(k,'=0 at M=',sp.solve(sp.Eq(s[k],0),M))
# M->oo limits
print({k:sp.limit(s[k],M,sp.oo) for k in [aA,aB,aP]})
# M=4m ceiling force
d={M:4*m}
F=sp.simplify(2*s[T2].subs(d)); print('ceiling F(M=4m)=',F)
acc={k:sp.simplify(s[k].subs(d)) for k in [aA,aB,aP]}; print(acc)
# COM check: F - 7mg = m aA + 2m aB + 4m aC, aC=-aP
assert sp.simplify(F-7*m*g-(m*acc[aA]+2*m*acc[aB]-4*m*acc[aP]))==0
# ceiling force < total weight for all M
expr=sp.simplify((M+3*m)*g-2*s[T2]); print('weight - ceiling =',sp.factor(expr))
# tensions positive
assert all(sp.simplify(s[T1].subs(M,v*m))>0 for v in [0.1,1,5,100])
print('OK')
