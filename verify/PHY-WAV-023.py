# PHY-WAV-023 検証：おもりを等間隔につけた弦（ビーズ弦）の固有振動
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp
S,m,b,k,t,w=sp.symbols('S m b k t omega',positive=True)
j=sp.symbols('j',integer=True)
A=sp.symbols('A')
y=lambda jj: A*sp.sin(jj*k*b)*sp.cos(w*t)
# (1) 運動方程式 m y_j'' = (S/b)(y_{j+1}-2y_j+y_{j-1})
lhs=m*sp.diff(y(j),t,2); rhs=S/b*(y(j+1)-2*y(j)+y(j-1))
res=sp.simplify(sp.expand_trig(sp.expand(lhs-rhs)).subs(w,2*sp.sqrt(S/(m*b))*sp.sin(k*b/2)))
print('(1) residual with dispersion:',sp.simplify(sp.trigsimp(res)))
# (2)(3) n=3 の行列固有値
for n in [2,3,4,5]:
    K=np.zeros((n,n))
    for i in range(n):
        K[i,i]=2
        if i>0:K[i,i-1]=-1
        if i<n-1:K[i,i+1]=-1
    ev,vec=np.linalg.eigh(K)   # omega^2 = (S/mb) ev
    om=np.sqrt(ev)
    formula=[2*np.sin(l*np.pi/(2*(n+1))) for l in range(1,n+1)]
    print('n=',n,'numeric',np.round(om,6),'formula',np.round(formula,6),' f2/f1=',om[1]/om[0])
# n=3 のモード形
n=3
K=np.array([[2,-1,0],[-1,2,-1],[0,-1,2]],float)
ev,vec=np.linalg.eigh(K)
print('mode shapes (columns):\n',np.round(vec/np.abs(vec).max(axis=0),4))
print('sqrt(2+sqrt2)=',np.sqrt(2+np.sqrt(2)))
# (4) 位相速度の比
v1=2*np.sin(np.pi/8)/(np.pi/(4)); v3=2*np.sin(3*np.pi/8)/(3*np.pi/4)
print('(4) v3/v1 =',v3/v1,' (1+sqrt2)/3 =',(1+np.sqrt(2))/3)
# 総質量・長さを同じにして小球の数を変えたときの基本振動数（一様な弦は π）
def f1(n): # n beads, total mass Mtot=1, L=1, S=1
    mm=1/n; bb=1/(n+1); return 2*np.sqrt(1/(mm*bb))*np.sin(np.pi/(2*(n+1)))
print('(5) n=3 f1',f1(3),' n=7',f1(7),' cont(mass1 uniform) pi',np.pi)
# (6) 新設問：R_n = (2/π) sqrt(n(n+1)) sin(π/(2(n+1)))，R_n<1，1-R_n ≈ 1/(2n)
nn_=sp.symbols('n',positive=True)
Rn=2/sp.pi*sp.sqrt(nn_*(nn_+1))*sp.sin(sp.pi/(2*(nn_+1)))
for nv in range(1,60):
    # 行列の固有値から独立に：小球 n 個，総質量 1，L=1，S=1
    K=2*np.eye(nv)-np.eye(nv,k=1)-np.eye(nv,k=-1)
    mm=1/nv; bb=1/(nv+1)
    w1=np.sqrt(np.linalg.eigvalsh(K)[0]/(mm*bb))
    R=w1/np.pi
    assert abs(R-float(Rn.subs(nn_,nv)))<1e-10 and R<1
print('R_3 =',float(Rn.subs(nn_,3)),' 4sqrt3/pi sin(pi/8)=',4*np.sqrt(3)/np.pi*np.sin(np.pi/8))
eps=sp.symbols('epsilon',positive=True)   # eps = 1/n
ser=sp.series(Rn.subs(nn_,1/eps),eps,0,3); print('R_n series in 1/n:',ser)
assert sp.simplify(ser.removeO().coeff(eps,1)+sp.Rational(1,2))==0
# 分散の効果のみ: sin x / x （x=π/(2(n+1)))は 1/n^2 の量；密度の効果 sqrt(n/(n+1)) が 1/(2n)
disp=sp.series((sp.sin(sp.pi/(2*(1/eps+1)))/(sp.pi/(2*(1/eps+1)))),eps,0,3); print('dispersion factor:',disp)
assert disp.removeO().coeff(eps,1)==0
dens=sp.series(sp.sqrt((1/eps)/(1/eps+1)),eps,0,2); print('density factor:',dens)
for nv in [10,100,1000]:
    print('n=',nv,' 1-R_n=',1-float(Rn.subs(nn_,nv)),' 1/(2n)=',1/(2*nv))
# 数値積分で n=3 の第1モードを確認
def rhs(tt,Y):
    yy=Y[:3]; vv=Y[3:]
    full=np.r_[0,yy,0]
    acc=(full[2:]-2*full[1:-1]+full[:-2])  # S/(m b)=1
    return np.r_[vv,acc]
Y0=np.r_[np.sin(np.arange(1,4)*np.pi/4),0,0,0]
r=solve_ivp(rhs,[0,30],Y0,rtol=1e-10,atol=1e-12,dense_output=True)
T1=2*np.pi/(2*np.sin(np.pi/8))
print('after one period (mode1):',np.round(r.sol(T1)[:3],6),' initial',np.round(Y0[:3],6))
