# S0121-1 内部振動しているばね連結体に当たった球のはね返り（見かけの反発係数が1を超える）
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,V,u,k=sp.symbols('m V u k',positive=True)
# (1) 静止した D に C が衝突
vA=V; vC=0; vG=V/2; e1=vG/V; Eint=sp.Rational(1,2)*m*V**2-sp.Rational(1,2)*2*m*vG**2
print('(1) e=',e1,' internal',sp.simplify(Eint))
# (2) 自然長で A が -u, B が +u のとき衝突（右向き正，C は +V）
# 等質量の弾性衝突で速度交換: C -> -u, A -> V
vG2=(V+u)/2; e2=sp.simplify((vG2-(-u))/V); print('(2) vG=',vG2,' e=',e2)
print('e>1 <=>',sp.solve(sp.Eq(e2,1),u))
Ein_before=m*u**2; Ein_after=sp.simplify(2*sp.Rational(1,2)*m*((V-u)/2)**2)
tot_b=sp.Rational(1,2)*m*V**2+m*u**2; tot_a=sp.Rational(1,2)*m*u**2+sp.Rational(1,2)*m*V**2+sp.Rational(1,2)*m*u**2
print('energy conserved',sp.simplify(tot_b-tot_a),' internal after',Ein_after)
# 数値シミュレーション：硬いばね（接触）で C と A の衝突を表し，再衝突の有無を確認
def sim(Vv,uv,kv=1.0,mv=1.0,kc=1e5,T=60):
    w=np.sqrt(2*kv/mv); l=1.0
    # 初期: A at 0 (自然長), B at l, A vel -uv, B vel +uv; C at -0.05 with Vv
    # Cを接触直前に置き，Aの振動位相を合わせるため，C-A間の距離0から開始
    def f(t,Y):
        xc,xa,xb,vc,va,vb=Y
        ov=xc-xa; Fc=kc*ov if ov>0 else 0.0
        Fs=kv*((xb-xa)-l)
        return [vc,va,vb,-Fc/mv,(Fc+Fs)/mv,-Fs/mv]
    s=solve_ivp(f,[0,T],[0.0,0.0,l,Vv,-uv,uv],rtol=1e-10,atol=1e-12,max_step=1e-3)
    xc,xa,xb,vc,va,vb=s.y
    # 衝突終了後の C の速度
    sep=np.where(xc-xa<0)[0]
    i0=sep[0]
    later=np.where(xc[i0:]-xa[i0:]>0)[0]
    return vc[-1],va[i0],vb[i0],len(later)
for Vv,uv in [(1.0,0.2),(1.0,0.5),(1.0,1.5),(2.0,0.1)]:
    vcf,va0,vb0,nre=sim(Vv,uv)
    print('V,u',Vv,uv,' C final',vcf,'(theory',-uv,') A,B after',va0,vb0,' re-contacts',nre)
    assert abs(vcf+uv)<5e-3 and nre==0
