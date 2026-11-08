# O20142-1 検証：鉛直磁場中の斜面レールでおもりを引き上げる導体棒（浮き上がりの制約と損失最小の角）
from sympy import *
import numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import minimize_scalar
M,m,g,E,B,l,r,th=symbols('M m g E B ell r theta',positive=True)
c,s=cos(th),sin(th)
Iinf=(M-m*s)*g/(B*l*c)
vinf=(E-Iinf*r)/(B*l*c)
lam=simplify(Iinf*r/E)
# エネルギー収支：E I = I^2 r + (M - m sin)g v
assert simplify(E*Iinf-Iinf**2*r-(M-m*s)*g*vinf)==0
# 棒を押さえた場合の最小
th1=solve(diff(lam,th),th); print('crit',th1)
dl=simplify(diff(lam,th)); print('dlam/dth numerator sign ~', factor(simplify(dl*c**2*B*l*E/(g*r))))
Ninf=simplify(m*g*c-Iinf*B*l*s); print('Ninf',simplify(Ninf*c))
assert simplify(Ninf - g*(m-M*s)/c)==0
s1=m/M
assert simplify(Ninf.subs(th,asin(s1)))==0
lam1=simplify(lam.subs(th,asin(s1))); print('lam1',lam1)
assert simplify(lam1-g*r*sqrt(M**2-m**2)/(E*B*l))==0
N0=m*g*c-E*B*l/r*s
# 数値：Nが単調増加（電流単調減少）を確認，制約付き最小
vals=dict(M=3.0,m=1.0,g=9.8,E=40.0,B=1.0,l=0.5,r=0.4)
Mn,mn,gn,En,Bn,ln,rn=[vals[k] for k in ['M','m','g','E','B','l','r']]
k=mn*gn*rn/(En*Bn*ln); thc=np.arctan(k)
print('k',k,'thc',thc,'sin thc',np.sin(thc),'m/M',mn/Mn, 'cond EBl>gr sqrt(M^2-m^2):',En*Bn*ln>gn*rn*np.sqrt(Mn**2-mn**2))
def run(t_):
    cc,ss=np.cos(t_),np.sin(t_)
    def f(t,y):
        v=y[0]  # 棒の斜面下向き速さ＝おもりの上昇速さ
        I=(En-v*Bn*ln*cc)/rn
        return [(I*Bn*ln*cc+mn*gn*ss-Mn*gn)/(Mn+mn)]
    sol=solve_ivp(f,[0,150],[0.0],max_step=1e-2,rtol=1e-10,atol=1e-12)
    v=sol.y[0]; I=(En-v*Bn*ln*cc)/rn
    N=mn*gn*cc-I*Bn*ln*ss
    return sol,I,N
sol,I,N=run(thc*0.999)
assert np.all(np.diff(I)<=1e-12) and N.min()>=-1e-9
print('N0 at thc*0.999',N[0],'Ninf',N[-1])
sol,I,N=run(thc*1.01); print('slightly steeper: N0',N[0]); assert N[0]<0
# 制約付き最小の数値確認
lamf=lambda t_: (Mn-mn*np.sin(t_))*gn*rn/(En*Bn*ln*np.cos(t_))
res=minimize_scalar(lamf,bounds=(0,thc),method='bounded',options={'xatol':1e-10})
lamc=k/mn*(Mn*np.sqrt(1+k**2)-mn*k)
print('numeric min at',res.x,'thc',thc,'lam',res.fun,'formula',lamc)
assert abs(res.x-thc)<1e-4 and abs(res.fun-lamc)<1e-8
# 定常での I,v の数値照合
cc=np.cos(thc*0.999); sol,I,N=run(thc*0.999)
Ith=float(Iinf.subs({M:Mn,m:mn,g:gn,B:Bn,l:ln,th:thc*0.999}))
assert abs(I[-1]-Ith)<1e-6
print('OK')
