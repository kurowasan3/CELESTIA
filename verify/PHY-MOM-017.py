# K20211-1 verification: repeated elastic collisions with upward-moving hitters
import sympy as sp
m,M,u,V0,v,w=sp.symbols('m M u V0 v w',positive=True)
# axis: up positive. A incoming downward speed u -> velocity -u; hitter +V0
va,vh=sp.symbols('va vh')
sol=sp.solve([sp.Eq(m*(-u)+M*V0,m*va+M*vh),sp.Eq(va-vh,-((-u)-V0))],[va,vh],dict=True)[0]
print('A after:',sp.simplify(sol[va]),' hitter after:',sp.simplify(sol[vh]))
# hitter stops: vh=0
Ms=sp.solve(sp.Eq(sol[vh],0),M)[0]; print('M for stop',sp.simplify(Ms),' A speed',sp.simplify(sol[va].subs(M,Ms)))
# sequence with stop scheme, u0=V0
uu=1; Ml=[]
for n in range(1,6):
    Mn=Ms.subs({u:uu,V0:1,m:1}); vn=sp.simplify(sol[va].subs(M,Ms).subs({u:uu,V0:1}))
    Ml.append((n,Mn,vn)); uu=vn
print(Ml)  # expect M=2n+1, v=n+1
# monotonicity in M and u
f=sol[va]
print('df/dM =',sp.factor(sp.diff(f,M)),' df/du=',sp.factor(sp.diff(f,u)))
# greedy with M=K m
K,n=sp.symbols('K n',positive=True)
r=(K-1)/(K+1)
vn=K-(K-1)*r**n
fK=sp.simplify(f.subs({M:K*m,V0:1}))
print('check recursion:',sp.simplify(fK.subs(u,vn)-vn.subs(n,n+1)), ' v0=',vn.subs(n,0))
# K=9 numbers
import math
for N in range(0,9):
    vv=9-8*0.8**N; print(N, round(vv,4), round(vv**2,3))
# brute-force: any strategy with M<=9m cannot beat greedy (random search)
import random
def step(u,Mr): return ((Mr-1)*u+2*Mr)/(Mr+1)
best=0
for t in range(20000):
    uu=1
    for k in range(4): uu=step(uu,random.uniform(0.1,9))
    best=max(best,uu)
print('random best after 4:',best,' greedy:',9-8*0.8**4)
