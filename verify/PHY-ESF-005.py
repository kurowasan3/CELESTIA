# P56-1 鏡像電荷：接地導体球と点電荷
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp
k,q,R,D,m,th,r=sp.symbols('k q R D m theta r',positive=True)
qp,b=sp.symbols('qp b',positive=True)
# (1) 球面上 V=0 となる像電荷 -qp at distance b
# 点 (R cos th, R sin th) と (D,0),(b,0) の距離
d1=sp.sqrt(R**2+D**2-2*R*D*sp.cos(th)); d2=sp.sqrt(R**2+b**2-2*R*b*sp.cos(th))
sol=sp.solve([sp.Eq(q/(D-R),qp/(R-b)), sp.Eq(q/(D+R),qp/(R+b))],[qp,b],dict=True)
print('(1)',sol)
s=sol[0]
V=k*q/d1-k*s[qp]/d2.subs(b,s[b])
print('V on sphere check:',[sp.N(V.subs({k:1,q:1,R:1,D:3,th:t})) for t in [0.3,1.2,2.5]])
# (2) 誘導電荷 = 像電荷（ガウス）: 数値で表面電荷を積分
# 表面電荷密度 sigma = eps0 E_r(外側) ; 4πk eps0=1 → sigma = E_r/(4πk)
import math
kk,qq,RR,DD=1.0,1.0,1.0,3.0
qi=qq*RR/DD; bi=RR**2/DD
def Er(t):
    x,y=RR*math.cos(t),RR*math.sin(t)
    tot=0
    for (Q,x0) in [(qq,DD),(-qi,bi)]:
        dx,dy=x-x0,y; dd=math.hypot(dx,dy)
        tot+=kk*Q*(dx*math.cos(t)+dy*math.sin(t))/dd**3
    return tot
from scipy.integrate import quad
Qind=quad(lambda t: Er(t)/(4*math.pi*kk)*2*math.pi*RR**2*math.sin(t),0,math.pi)[0]
print('(2) induced charge numeric',Qind,' expected',-qq*RR/DD)
F=k*q*s[qp]/(D-s[b])**2
print('(2) F=',sp.simplify(F))
# (3) 無限遠から D まで運ぶ外力の仕事 W = -∫_inf^D F dD (引力なので外力は外向き F の逆)
Dp=sp.symbols('Dp',positive=True)
Fa=sp.integrate(F.subs(D,Dp),Dp)
W=sp.simplify(Fa.subs(Dp,D)-sp.limit(Fa,Dp,sp.oo))  # W_ext=∫_inf^D (+F)dD'
print('(3) W_ext from inf=',sp.simplify(W))
naive=-k*q*s[qp]/(D-s[b])
print('    naive q*V_image=',sp.simplify(naive),' ratio',sp.simplify(W/naive))
# 別解：エネルギー = (1/2) q V_image(位置) を確認（同じ式）
# (4) 3R -> 2R speed
U=sp.simplify(W)
v=sp.sqrt(2*(U.subs(D,3*R)-U.subs(D,2*R))/m)
print('(4) v=',sp.simplify(v))
# 数値積分
kk,qq,RR,mm=1.0,1.0,1.0,1.0
def f(t,y):
    x,vx=y
    return [vx,-kk*qq**2*RR*x/(x**2-RR**2)**2/mm]
ev=lambda t,y:y[0]-2*RR; ev.terminal=True
so=solve_ivp(f,[0,100],[3*RR,0],events=ev,rtol=1e-11,atol=1e-12)
print('numeric v at 2R',abs(so.y_events[0][0][1]),' formula',math.sqrt(5/24))
