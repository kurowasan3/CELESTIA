# J017-1 検証：液体の底にちょうつがいで留めた軽い棒の傾き（鉛直が安定か，傾くかの分岐）とちょうつがいの力
import sympy as sp, numpy as np
l,S,rho,rho0,g,d,th=sp.symbols('l S rho rho0 g d theta',positive=True)
# 棒が傾いて上端が液面上に出ているとき：液中の長さ s=d/sin(th)
s=d/sp.sin(th)
# ちょうつがい O まわりモーメント（反時計回りを正，棒は O から右上へ）：浮力は液中部分の中点 s/2，重力は l/2
tau=rho0*S*s*g*(s/2)*sp.cos(th)-rho*S*l*g*(l/2)*sp.cos(th)
sol=sp.solve(sp.Eq(rho0*s**2,rho*l**2),th); print('sin th solutions:',[sp.simplify(sp.sin(x)) for x in sol])
s_eq=l*sp.sqrt(rho/rho0); print('submerged length s =',s_eq)
# 鉛直状態（0<d<l）の安定性：わずかに傾けた（角度 phi = 90deg - th 小）ときのモーメントの向き
phi=sp.symbols('phi',positive=True)
s_v=d/sp.cos(phi)  # 鉛直からphi傾けたときの液中の長さ
restoring=(rho0*S*s_v*g*(s_v/2)*sp.sin(phi)-rho*S*l*g*(l/2)*sp.sin(phi))  # 浮力は鉛直に戻す向き・重力は倒す向き：正なら鉛直に戻す
lin=sp.series(restoring,phi,0,2).removeO(); print('restoring torque (linear):',sp.factor(lin))
print(' => vertical stable iff rho0*d^2 > rho*l^2')
# ちょうつがいの力（棒に及ぼす力，鉛直下向き正）
F_tilt=sp.simplify(rho0*S*s_eq*g-rho*S*l*g); print('hinge force tilted:',sp.factor(F_tilt))
F_vert=rho0*S*d*g-rho*S*l*g; F_sub=(rho0-rho)*S*l*g
dc=l*sp.sqrt(rho/rho0)
print('continuity at d=dc:',sp.simplify(F_vert.subs(d,dc)-F_tilt))
# 水平な力はない？（浮力・重力とも鉛直 → ちょうつがいの力も鉛直）
# 数値：ポテンシャルエネルギー最小で平衡角を確認
ln,rn,r0n,gn,Sn=1.0,0.36,1.0,9.8,1.0
def U(thv,dn):
    # 重力の位置エネルギー + 浮力の位置エネルギー（-ρ0 g × 液中体積の重心高さ...）
    # 浮力の位置エネルギー = +rho0 g ∫(液中) (液面までの深さ) dV（沈めるのに仕事が要る）
    sub=min(ln,dn/np.sin(thv))
    return rn*Sn*ln*gn*(ln/2)*np.sin(thv) + r0n*gn*Sn*(dn*sub - sub**2*np.sin(thv)/2)
for dn in [0.3,0.5,0.59,0.7,0.9]:
    ths=np.linspace(0.01,np.pi/2,200000); Us=[U(t_,dn) for t_ in ths[::100]]
    i=int(np.argmin(Us)); tmin=ths[::100][i]
    pred=np.degrees(np.arcsin(min(1,dn/(ln*np.sqrt(rn/r0n)))))
    print('d=',dn,' numeric eq angle',np.degrees(tmin),' predicted',pred)

# --- (4) 傾いた状態の安定性と，d_c 直下での傾き phi ≈ sqrt(2(d_c-d)/d_c)
Mth=S*g*sp.cos(th)/2*(rho0*(d/sp.sin(th))**2-rho*l**2)     # 起こす向きのモーメント
th_eq=sp.asin(d/dc)
dM=sp.simplify(sp.diff(Mth,th).subs(th,th_eq)); print('dM/dtheta at eq =',dM)
num={l:1,S:1,g:9.8,rho:0.36,rho0:1.0}
for dn in [0.1,0.3,0.5,0.59]:
    v_=float(dM.subs(num).subs(d,dn)); print(' d',dn,' dM/dtheta',v_); assert v_<0    # 復元的（安定）
    te=float(th_eq.subs(num).subs(d,dn))
    assert float(Mth.subs(num).subs(d,dn).subs(th,te+1e-3))<0 and float(Mth.subs(num).subs(d,dn).subs(th,te-1e-3))>0
dcn=0.6
for eps in [1e-2,1e-3,1e-4]:
    dn=dcn*(1-eps); phi_exact=np.pi/2-np.arcsin(dn/dcn); phi_app=np.sqrt(2*(dcn-dn)/dcn)
    print(' eps',eps,' phi exact',phi_exact,' approx',phi_app,' rel err',abs(phi_exact/phi_app-1)); assert abs(phi_exact/phi_app-1)<eps
print(' 2% below d_c: phi =',np.pi/2-np.arcsin(0.98),'rad =',np.degrees(np.pi/2-np.arcsin(0.98)),'deg')
print('(4) OK')
