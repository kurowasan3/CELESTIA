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
# ---- (4)(a) total hitter mass lower bound N(N+2)m ; stop scheme attains it ----
N=sp.symbols('N',positive=True,integer=True)
kk=sp.symbols('k',positive=True,integer=True)
print('(4a) sum (2k+1) k=1..N =',sp.factor(sp.summation(2*kk+1,(kk,1,N))))
assert sp.simplify(sp.summation(2*kk+1,(kk,1,N))-N*(N+2))==0
# energy bookkeeping per collision: KE_A gain = KE_hitter loss (elastic)
gainA=sp.simplify((m*sol[va]**2/2-m*u**2/2)-(M*V0**2/2-M*sol[vh]**2/2)); print('(4a) energy check per collision:',gainA)
assert gainA==0
# random strategies: (sum M)/m >= (v_final/V0)^2-1
import random
def step2(u_,Mr):  # m=1,V0=1
    vA=((Mr-1)*u_+2*Mr)/(Mr+1); vH=((Mr-1)-2*u_)/(Mr+1); return vA,vH
worst=1e9
for t_ in range(20000):
    uu=1.0; Ms_=0.0
    for k_ in range(random.randint(1,6)):
        Mr=random.uniform(0.1,30); uu,_=step2(uu,Mr); Ms_+=Mr
    worst=min(worst,Ms_-(uu**2-1))
print('(4a) min over random strategies of sumM-(v^2-1) =',worst); assert worst>=-1e-9
# ---- (4)(b) single collision bound 3V0 ----
v1=sp.simplify(sol[va].subs(u,V0)); print('(4b) v1 =',v1,' limit M->oo:',sp.limit(v1,M,sp.oo))
assert sp.simplify(v1-(3*V0-4*m*V0/(M+m)))==0
# (3) check: with M<=Km (K=9), no strategy beats greedy M=Km (random search)
def step(u,Mr): return ((Mr-1)*u+2*Mr)/(Mr+1)
best=0
for t_ in range(20000):
    uu=1
    for k_ in range(4): uu=step(uu,random.uniform(0.1,9))
    best=max(best,uu)
print('(3) random best after 4:',best,' greedy:',9-8*0.8**4); assert best<=9-8*0.8**4+1e-12
