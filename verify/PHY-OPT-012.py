# P39-1 液面の上昇に伴う見かけの大きさ（視角）の変化
import sympy as sp, numpy as np
from scipy.optimize import brentq
L,H,D,n=sp.symbols('L H D n',positive=True)
th_below=L/(H-D+D/n); th_above=n*L/H
print('continuity at D=H:',sp.simplify(th_below.subs(D,H)-th_above))
# 数値データからの逆算
inv0,inv1,invinf,D1=40.0,32.5,30.0,0.30
nn=inv0/invinf; LL=D1*(1-1/nn)/(inv0-inv1); HH=inv0*LL
print('n=',nn,' L=',LL,' H=',HH)
# 厳密な光線追跡で近軸近似を確認：端Aから出て眼に入る光線の角度
def exact_theta(Lv,Hv,Dv,nv):
    if Dv<Hv:
        f=lambda a: Dv*np.tan(a)+(Hv-Dv)*np.tan(np.arcsin(nv*np.sin(a)))-Lv
        a=brentq(f,0,np.arcsin(1/nv)-1e-9); return np.arcsin(nv*np.sin(a))
    else:
        return np.arctan(Lv/Hv)*1.0  # 眼の窓の直前まで液体，窓で屈折
for Dv in [0,0.1,0.3,0.39]:
    ex=exact_theta(LL,HH,Dv,nn); ap=LL/(HH-Dv+Dv/nn)
    print(Dv,1/ex,1/ap)
# D>H: 窓（空気側）での角 = asin(n sin(atan(L/H)))
a=np.arctan(LL/HH); print('above:',1/np.arcsin(nn*np.sin(a)), HH/(nn*LL))
# (4) 眼の高さ h だけ変えた場合：L/θ = h + D/n -> 2測定の差は h の差だけ
h1,h2=sp.symbols('h1 h2',positive=True)
t1=L/(h1+D/n); t2=L/(h2+D/n)
print(sp.simplify(L/t1-L/t2))
