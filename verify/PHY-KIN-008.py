# J002-1 検証：反応の遅れとブレーキ性能の違う2台の車の最小車間距離（場合分け）
import sympy as sp, numpy as np
V,tau,aA,aB,t=sp.symbols('V tau a_A a_B t',positive=True)
# 一般式
Dfinal=V*tau+V**2/(2*aB)-V**2/(2*aA)
Dmid=aA*aB*tau**2/(2*(aB-aA))
tstar=aB*tau/(aB-aA)
aBc=sp.solve(sp.Eq(tstar,V/aA),aB); print('boundary aB =',aBc)
# 境界で2式一致するか
print('Dmid-Dfinal at boundary:',sp.simplify((Dmid-Dfinal).subs(aB,aBc[0])))
num={V:25,tau:1,aA:6}
print('boundary numeric',sp.N(aBc[0].subs(num)))
print('aB=5 Dfinal',sp.N(Dfinal.subs(num).subs(aB,5)))
print('aB=6 Dfinal',sp.N(Dfinal.subs(num).subs(aB,6)))
print('aB=9 Dmid',sp.N(Dmid.subs(num).subs(aB,9)),' Dfinal',sp.N(Dfinal.subs(num).subs(aB,9)))
# 数値：位置を時刻で細かく計算して最小の (xA - xB) を求める（初期車間 0 とし，必要な D = -min）
def need(Vn,taun,aAn,aBn):
    ts=np.linspace(0,30,3000001)
    vA=np.maximum(Vn-aAn*ts,0); 
    xA=np.where(ts<Vn/aAn, Vn*ts-0.5*aAn*ts**2, Vn**2/(2*aAn))
    tb=np.maximum(ts-taun,0)
    xB=np.where(ts<taun, Vn*ts, Vn*taun+np.where(tb<Vn/aBn, Vn*tb-0.5*aBn*tb**2, Vn**2/(2*aBn)))
    gap=xA-xB; i=gap.argmin(); return -gap[i], ts[i]
for aBn in [5,6,7.5,7.9,8.5,9,12]:
    print('aB=',aBn,' numeric need, t_at =',need(25,1,6,aBn))
# V 依存性：aB=9 で V を変える
for Vn in [20,25,30]:
    print('V=',Vn,' aB=9 need',need(Vn,1,6,9), ' boundary aB',Vn*6/(Vn-6))
