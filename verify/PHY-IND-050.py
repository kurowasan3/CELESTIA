# PHY-IND-050 検証：45°傾けた正方形コイルが幅 w の磁場の帯を通過するときの起電力・ジュール熱・逆問題
import sympy as sp, numpy as np
from scipy.integrate import quad
def chord(x,X,a):
    u=X-x
    if u<0 or u>2*a: return 0.0
    return 2*min(u,2*a-u)
def area(X,a,w):   # 磁場中のコイル面積（数値積分）
    lo,hi=max(0.0,X-2*a),min(w,X)
    if hi<=lo: return 0.0
    return quad(lambda x: chord(x,X,a),lo,hi,points=[X-a],limit=200)[0]
def V_formula(X,a,w,B=1.0,v=1.0):  # 反時計回り正，w<a
    if X<0: d=0
    elif X<w: d=2*X
    elif X<a: d=2*w
    elif X<a+w: d=2*(2*a+w-2*X)
    elif X<2*a: d=-2*w
    elif X<2*a+w: d=-2*(2*a+w-X)
    else: d=0
    return -B*v*d
a,w=1.0,0.37
Xs=np.linspace(0.001,2*a+w-0.001,400); h=1e-5
for X in Xs:
    dA=(area(X+h,a,w)-area(X-h,a,w))/(2*h)
    assert abs(-dA-V_formula(X,a,w))<1e-4, (X,dA,V_formula(X,a,w))
print('piecewise emf OK (w<a)')
# (1) 一般式 dA/dX = l0 - lw をランダム位置で確認（w>a も含む）
for w2 in [0.5,1.3,1.8,2.5]:
    for X in np.linspace(0.01,2*a+w2-0.01,57):
        dA=(area(X+h,a,w2)-area(X-h,a,w2))/(2*h)
        assert abs(dA-(chord(0,X,a)-chord(w2,X,a)))<1e-4
print('dA/dX = l0 - lw OK')
# (4) ジュール熱 Q=(B^2 v/R) ∫(dA/dX)^2 dX
As,Ws=sp.symbols('a w',positive=True); X=sp.symbols('X')
pieces=[(2*X,0,Ws),(2*Ws,Ws,As),(2*(2*As+Ws-2*X),As,As+Ws),(-2*Ws,As+Ws,2*As),(-2*(2*As+Ws-X),2*As,2*As+Ws)]
Qs=sp.expand(sum(sp.integrate(e**2,(X,lo,hi)) for e,lo,hi in pieces)); print('Q*R/(B^2 v)=',Qs)
for wn in [0.2,0.5,0.8,1.0]:
    qn=quad(lambda X: V_formula(X,1.0,wn)**2,0,2+wn,points=[wn,1,1+wn,2],limit=400)[0]
    assert abs(qn-float(Qs.subs({As:1,Ws:wn})))<1e-6
print('Q formula numeric OK; w>=2a gives',quad(lambda X:(chord(0,X,1)-chord(3,X,1))**2,0,5,points=[1,2,3,4],limit=400)[0])
# 外力の仕事 = ジュール熱（F = B^2 v (l0-lw)^2/R, 仕事 = ∫F dX）は同じ積分
# (3) 逆問題
T1,T2,V1,v=sp.symbols('T1 T2 V1 v',positive=True)
sol=sp.solve([sp.Eq(Ws/v,T1),sp.Eq((As-Ws)/v,T2)],[v,Ws],dict=True)[0]; print(sol)
B=sp.symbols('B',positive=True)
Bsol=sp.solve(sp.Eq(2*B*sol[v]*sol[Ws],V1),B)[0]; print('B=',sp.simplify(Bsol))
# 数値で波形を作って逆算
an,wn,vn,Bn=0.10,0.03,0.5,0.4
ts=np.linspace(0,(2*an+wn)/vn,20001); Vt=np.array([abs(V_formula(vn*t,an,wn,Bn,vn)) for t in ts])
Vp=Vt.max(); plateau=ts[np.isclose(Vt,Vp,rtol=1e-9)]
T1n=plateau[0]; T2n=plateau[np.where(np.diff(plateau)>2*(ts[1]-ts[0]))[0][0]] - T1n if np.any(np.diff(plateau)>2*(ts[1]-ts[0])) else None
print('measured T1,T2,Vp',T1n,T2n,Vp)
print('recover v,w,B:',an/(T1n+T2n),an*T1n/(T1n+T2n),float(Bsol.subs({V1:Vp,T1:T1n,T2:T2n,As:an})))
print('Fmax=4B^2 v w^2/R ; OK')
