# P18 類題：楕円軌道の近地点と遠地点での噴射（オーベルト効果）
import sympy as sp, numpy as np
g,R=sp.symbols('g R',positive=True)
GM=g*R**2; rp=2*R; ra=6*R
vp,va=sp.symbols('vp va',positive=True)
s=sp.solve([sp.Eq(vp*rp,va*ra),sp.Eq(vp**2/2-GM/rp,va**2/2-GM/ra)],[vp,va],dict=True)[0]
print('vp=',sp.simplify(s[vp]),' va=',sp.simplify(s[va]))
dvp=sp.sqrt(2*GM/rp)-s[vp]; dva=sp.sqrt(2*GM/ra)-s[va]
print('dvp=',sp.simplify(dvp),' dva=',sp.simplify(dva),' ratio',sp.nsimplify(sp.simplify(dva/dvp)),float(dva/dvp))
sq=np.sqrt(9.8*6.4e6)/1000
print('sqrt(gR)=',sq,' dvp=',float(dvp/sp.sqrt(g*R))*sq,' dva=',float(dva/sp.sqrt(g*R))*sq)
# 同じ dv での比エネルギー増加
dv=sp.symbols('dv',positive=True)
dEp=sp.expand((s[vp]+dv)**2/2-s[vp]**2/2); dEa=sp.expand((s[va]+dv)**2/2-s[va]**2/2)
print('dE perigee',dEp,' apogee',dEa)
# ロケット＋ガスのエネルギー
M,dm,v,w=sp.symbols('M dm v w',positive=True)
Dv=dm*w/M
K=sp.Rational(1,2)*(M-dm)*(v+Dv)**2+sp.Rational(1,2)*dm*(v+Dv-w)**2-sp.Rational(1,2)*M*v**2
print('total dK =',sp.simplify(K))
assert sp.diff(sp.simplify(K),v)==0
Krocket=sp.Rational(1,2)*(M-dm)*((v+Dv)**2-v**2)
Kgas=sp.Rational(1,2)*dm*((v+Dv-w)**2-v**2)
print('rocket gain',sp.expand(Krocket),'\n gas change',sp.expand(Kgas))
# 数値軌道で：近地点で dvp だけ加速→エネルギー0（放物線）
gn,Rn=9.8,6.4e6; GMn=gn*Rn**2
vpn=float(s[vp].subs({g:gn,R:Rn})); dvpn=float(dvp.subs({g:gn,R:Rn}))
E=(vpn+dvpn)**2/2-GMn/(2*Rn); print('energy after perigee burn',E)
van=float(s[va].subs({g:gn,R:Rn})); dvan=float(dva.subs({g:gn,R:Rn}))
print('energy after apogee burn',(van+dvan)**2/2-GMn/(6*Rn))
