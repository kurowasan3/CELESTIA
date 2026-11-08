# P57-1 トムソン模型：一様帯電球内の2電子・3電子
import sympy as sp, numpy as np, math
from scipy.integrate import solve_ivp
k,e,R,m,a,r,s=sp.symbols('k e R m a r s',positive=True)
# (1) ガウス: E*4πr^2 = 4πk * 2e r^3/R^3
E=sp.solve(sp.Eq(sp.Symbol('E')*4*sp.pi*r**2,4*sp.pi*k*2*e*r**3/R**3),sp.Symbol('E'))[0]
print('(1)E=',E)
# (2) 2電子のつり合い
a2=sp.solve(sp.Eq(e*E.subs(r,a),k*e**2/(2*a)**2),a)
print('(2) a=',a2)
# (3) x軸上の力：電子1(位置x1),電子2(位置x2) ; F1 = -2ke^2 x1/R^3 + ke^2/(x1-x2)^2 (x1>x2)
x1,x2=sp.symbols('x1 x2',real=True)
F1=-2*k*e**2*x1/R**3+k*e**2/(x1-x2)**2
F2=-2*k*e**2*x2/R**3-k*e**2/(x1-x2)**2
J=sp.Matrix([[sp.diff(F1,x1),sp.diff(F1,x2)],[sp.diff(F2,x1),sp.diff(F2,x2)]]).subs({x1:R/2,x2:-R/2})
ev=(-J/m).eigenvals()
print('(3) omega^2:',[sp.simplify(v) for v in ev])
# 数値：呼吸モードの周期
kk,ee,RR,mm=1,1,1,1
def f(t,y):
    X1,X2,V1,V2=y
    f1=-2*X1+1/(X1-X2)**2; f2=-2*X2-1/(X1-X2)**2
    return [V1,V2,f1,f2]
eps=1e-3
so=solve_ivp(f,[0,20],[0.5+eps,-0.5-eps,0,0],rtol=1e-12,atol=1e-14,dense_output=True)
t=np.linspace(0,20,200001); X=so.sol(t)[0]-0.5
zc=t[1:][(X[:-1]>0)&(X[1:]<=0)]
print('breathing period numeric',zc[1]-zc[0],' theory',2*math.pi/math.sqrt(6))
so=solve_ivp(f,[0,20],[0.5+eps,-0.5+eps,0,0],rtol=1e-12,atol=1e-14,dense_output=True)
X=so.sol(t)[0]-0.5; zc=t[1:][(X[:-1]>0)&(X[1:]<=0)]
print('COM period numeric',zc[1]-zc[0],' theory',2*math.pi/math.sqrt(2))
# (4) 3電子 (球の電荷3e) 正三角形
E3=k*3*e*a/R**3
rep=2*k*e**2/(3*a**2)*sp.cos(sp.pi/6)
a3=sp.solve(sp.Eq(e*E3,rep),a)
print('(4) a=',a3,' side=',[sp.simplify(sp.sqrt(3)*x) for x in a3])
# 数値で3電子の平衡を緩和法で確認
def force(P):
    F=-3*P.copy()
    for i in range(3):
        for j in range(3):
            if i!=j:
                d=P[i]-P[j]; F[i]+=d/np.linalg.norm(d)**3
    return F
P=np.array([[0.3,0,0.01],[-0.2,0.4,0],[-0.1,-0.3,-0.02]])
for _ in range(20000): P+=0.01*force(P)
print('numeric distances',[np.linalg.norm(P[i]-P[j]) for i,j in [(0,1),(1,2),(0,2)]],' from center',np.linalg.norm(P,axis=1))
