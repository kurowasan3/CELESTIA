# P13-1: 台のはり出した鉛直ループを回る小球と台の転倒
import sympy as sp, numpy as np
th,lam,m,g,R=sp.symbols('theta lambda m g R',positive=True)
c=sp.cos(th); s=sp.sin(th)
N=m*g*(lam-2+3*c)
# ベクトル外積で C まわりのモーメント（時計回りを正）を独立に計算
O=sp.Matrix([0,R]); C=sp.Matrix([0,0]); P=O+R*sp.Matrix([s,-c])
F=N*sp.Matrix([s,-c])           # 小球が軌道を押す力（O→小球の向き）
r=P-C
tau_ccw=sp.simplify(r[0]*F[1]-r[1]*F[0]); print('torque cw =',sp.simplify(-tau_ccw))
Mom=-tau_ccw
d=sp.diff(Mom.subs(lam,5),th); print('stationary lam=5:',sp.solve(sp.Eq(sp.simplify(d),0),th))
print('max lam=5:',sp.nsimplify(sp.simplify(Mom.subs({lam:5,th:sp.pi/3})/(m*g*R))))
# lam=9: cos=1/3
val=sp.simplify(Mom.subs(lam,9).subs(th,sp.acos(sp.Rational(1,3)))/(m*g*R)); print('lam=9 at cos=1/3:',val)
# 数値で max を直接探索
for L_ in [5,9]:
    f=sp.lambdify(th,(Mom/(m*g*R)).subs(lam,L_)); x=np.linspace(0,np.pi,200001); y=f(x)
    print(L_, 'num max',y.max(),'at cos',np.cos(x[y.argmax()]))
# 一般の lam: cos th* = (-(lam-2)+sqrt((lam-2)^2+72))/12
for L_ in [5,6,9,20,100]:
    cs=(-(L_-2)+np.sqrt((L_-2)**2+72))/12
    f=sp.lambdify(th,(Mom/(m*g*R)).subs(lam,L_)); x=np.linspace(0,np.pi,200001)
    print(L_,cs,np.cos(x[f(x).argmax()]))
# 最大モーメントが lam について単調増加か
Fm=[]
for L_ in np.linspace(5,12,15):
    f=sp.lambdify(th,(Mom/(m*g*R)).subs(lam,L_)); Fm.append(f(np.linspace(0,np.pi,20001)).max())
print('monotone:',all(np.diff(Fm)>0))
