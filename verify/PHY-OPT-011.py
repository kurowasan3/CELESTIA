# P38-1 主虹・副虹の角度，停留条件，色の順，存在条件
import numpy as np, sympy as sp
from scipy.optimize import minimize_scalar
a,n=sp.symbols('alpha n',positive=True)
b=sp.asin(sp.sin(a)/n)
Th1=4*b-2*a                 # 対日点方向からの角（主虹）
Th2=sp.pi+2*a-6*b           # 副虹
for k,Th in [(1,Th1),(2,Th2)]:
    d=sp.diff(Th,a)
    # 停留条件の確認：cos^2 a = (n^2-1)/(k(k+2))
    c2=(n**2-1)/(k*(k+2))
    a0=sp.acos(sp.sqrt(c2))
    for nv in [sp.Rational(4,3), sp.Rational(3,2)]:
        print(k, nv, sp.N(d.subs(a,a0).subs(n,nv)))
# 直接の光線追跡（ベクトル）で偏角を計算し，式と照合
def trace(bimp,nn,k):
    # 単位円，光は +x 向き，衝突径数 bimp (y座標 = -bimp で下から，符号は無関係)
    p=np.array([-np.sqrt(1-bimp**2),bimp]); d=np.array([1.,0.])
    def refract(d,nrm,n1,n2):
        cosi=-np.dot(nrm,d); r=n1/n2; k2=1-r*r*(1-cosi**2)
        return r*d+(r*cosi-np.sqrt(k2))*nrm
    nrm=p/np.linalg.norm(p); d=refract(d,nrm,1,nn)
    for i in range(k+1):
        # 球内を進む
        t=-2*np.dot(p,d); p=p+t*d
        nrm=p
        if i<k: d=d-2*np.dot(d,nrm)*nrm
    d=refract(d,-nrm,nn,1) if False else None
    # 出射
    nrm_out=p; cosi=np.dot(nrm_out,dprev) if False else None
    return p
def deviation(bimp,nn,k):
    p=np.array([-np.sqrt(1-bimp**2),bimp]); d=np.array([1.,0.])
    def refr(d,nrm,r):
        cosi=-np.dot(nrm,d); k2=1-r*r*(1-cosi**2)
        return r*d+(r*cosi-np.sqrt(k2))*nrm
    d=refr(d,p,1/nn)
    for i in range(k+1):
        t=-2*np.dot(p,d); p=p+t*d
        if i<k: d=d-2*np.dot(d,p)*p
    d=refr(d,-p,nn)
    # 対日点方向（-x）とのなす角
    return np.degrees(np.arccos(np.clip(-d[0],-1,1)))
def formula(bimp,nn,k):
    al=np.arcsin(bimp); be=np.arcsin(bimp/nn)
    return np.degrees(4*be-2*al) if k==1 else np.degrees(abs(np.pi+2*al-6*be))
for bb in [0.3,0.7,0.9]:
    print('trace vs formula', deviation(bb,4/3,1),formula(bb,4/3,1), deviation(bb,4/3,2),formula(bb,4/3,2))
for nn in [4/3,1.331,1.343,1.5]:
    r1=minimize_scalar(lambda x:-deviation(x,nn,1),bounds=(0,0.9999),method='bounded')
    r2=minimize_scalar(lambda x:deviation(x,nn,2),bounds=(0,0.9999),method='bounded')
    c1=np.sqrt((nn*nn-1)/3); c2=np.sqrt((nn*nn-1)/8)
    print(nn,'primary max %.3f'%(-r1.fun),'formula %.3f'%formula(np.sin(np.arccos(c1)),nn,1),
          ' secondary min %.3f'%r2.fun,'formula %.3f'%formula(np.sin(np.arccos(c2)),nn,2))
print('alpha,beta n=4/3 k=2', np.degrees(np.arccos(np.sqrt((16/9-1)/8))), np.degrees(np.arcsin(np.sin(np.arccos(np.sqrt((16/9-1)/8)))*3/4)))
print('alpha,beta n=4/3 k=1', np.degrees(np.arccos(np.sqrt((16/9-1)/3))), np.degrees(np.arcsin(np.sin(np.arccos(np.sqrt((16/9-1)/3)))*3/4)))
# n=2.2 で主虹（極大）がないこと：Th1 は単調増加
bs=np.linspace(0.001,0.999,2000); v=[formula(x,2.2,1) for x in bs]
print('n=2.2 monotonic:',np.all(np.diff(v)>0))
v=[formula(x,2.5,2) for x in bs]; print('n=2.5 secondary has interior min:', np.argmin(v) not in (0,len(v)-1))
v=[formula(x,3.2,2) for x in bs]; print('n=3.2 secondary min at end?:', np.argmin(v), len(v))
