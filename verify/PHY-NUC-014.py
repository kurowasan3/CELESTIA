# P89-1: 静止したポジトロニウムの3光子消滅の運動学
import numpy as np
mc2=511.0
# (2) 正弦定理: E_k ∝ sin φ_ij
phi12,phi23,phi31=np.radians([90,135,135])
w=np.array([np.sin(phi23),np.sin(phi31),np.sin(phi12)]); E=2*mc2*w/w.sum(); print('E1,E2,E3',E, E/mc2)
# 運動量和の確認: 方向 θ1=0, θ2=90°, θ3=90+135=225°
th=np.radians([0,90,225]); p=E[:,None]*np.c_[np.cos(th),np.sin(th)]; print('sum p',p.sum(0))
# (3) 各光子 E<=mc^2 : ランダムな運動学的に許される配置を生成して最大値を確認
rng=np.random.default_rng(1); mx=0
for _ in range(200000):
    a=rng.uniform(0,2*np.pi,3)
    d=np.c_[np.cos(a),np.sin(a)]
    # E1 d1+E2 d2+E3 d3=0 の非負解（E1+E2+E3=2mc2）を求める
    A=np.vstack([d.T,np.ones(3)]); b=np.array([0,0,2*mc2])
    try: x=np.linalg.solve(A,b)
    except: continue
    if (x>=0).all(): mx=max(mx,x.max())
print('max single photon energy found',mx)
# (4) E1=0.80 mc^2 のとき，光子2,3のなす角の最小値
E1=0.8*mc2; S=2*mc2-E1
E2=np.linspace(E1/2+ (S-E1)/2*0+1e-9, S-1e-9,200001)  # 条件 |E2-E3|<=E1 の範囲
E2=E2[np.abs(2*E2-S)<=E1]; E3=S-E2
cospsi=(E1**2-E2**2-E3**2)/(2*E2*E3)
psi=np.degrees(np.arccos(np.clip(cospsi,-1,1)))
print('psi min',psi.min(),' at E2=',E2[np.argmin(psi)],' formula',2*np.degrees(np.arccos(E1/S)))
print('psi max',psi.max())
