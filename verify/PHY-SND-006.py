# P37-1 風のある上層での音の屈折：分散関係（位相整合）による独立導出と数値
import sympy as sp
import numpy as np
v1,v2,w,th,ph,om=sp.symbols('v1 v2 w theta phi omega',positive=True)
kx=om*sp.sin(th)/v1                      # 境界に沿う波数は保存
# 上層：媒質の速度 (w,0) ，分散関係 omega - w*kx = v2*|k|
# 風下: 上層の波面法線 (sin ph, cos ph), |k| = kx/sin(ph)
eq_down=sp.Eq(om - w*kx, v2*kx/sp.sin(ph))
sol=sp.solve(eq_down, sp.sin(ph))[0]
print('down sin(phi)=',sp.simplify(sol))
target=v2*sp.sin(th)/(v1-w*sp.sin(th))
print('down check:',sp.simplify(sol-target))
# Huygens 形式 v1/sin th = v2/sin ph + w
print('huygens down:',sp.simplify((v1/sp.sin(th)-w) - v2/sol))
# 風上: kx<0 ->  omega + w*|kx| = v2*|k|
eq_up=sp.Eq(om + w*kx, v2*kx/sp.sin(ph))
solu=sp.solve(eq_up, sp.sin(ph))[0]
print('up sin(phi)=',sp.simplify(solu), ' check:',sp.simplify(solu - v2*sp.sin(th)/(v1+w*sp.sin(th))))
# 群速度 = v2 * khat + w xhat を数値微分で確認
def omega(kxv,kzv,v2v,wv): return v2v*np.hypot(kxv,kzv)+wv*kxv
v2n,wn=330.,20.
k=np.array([0.3,0.8]); h=1e-7
g=[(omega(k[0]+h,k[1],v2n,wn)-omega(k[0]-h,k[1],v2n,wn))/(2*h),(omega(k[0],k[1]+h,v2n,wn)-omega(k[0],k[1]-h,v2n,wn))/(2*h)]
print('group vel',g,' expected',v2n*k/np.linalg.norm(k)+np.array([wn,0]))
# 風上で光線が鉛直になる条件: ray x成分 = -v2 sin ph + w = 0
sth=sp.symbols('s',positive=True)
s_up=sp.solve(sp.Eq(v2*sth/(v1+w*sth), w/v2), sth)[0]
print('upwind vertical ray sin(theta)=',sp.simplify(s_up))
# (4) 新設問：「上層は音速 v2+w の媒質」という誤った主張との比較
s_=sp.symbols('s',positive=True)   # s = sin(theta)
true_sphi=v2*s_/(v1-w*s_)
claim_sphi=(v2+w)*s_/v1
sols=sp.solve(sp.Eq(true_sphi,claim_sphi),s_)
print('(4) claim agrees only at sin(theta)=',sols)
assert len(sols)==1 and sp.simplify(sols[0]-v1/(v2+w))==0
assert sp.simplify(true_sphi.subs(s_,v1/(v2+w))-1)==0     # 臨界角で sin(phi)=1
# 法線方向の速さ v2 + w sin(phi) を使った交線速度の一致が (1) と同じ式を与える
sp_=sp.symbols('sp',positive=True)
sol_n=sp.solve(sp.Eq(v1/s_,(v2+w*sp_)/sp_),sp_)[0]
assert sp.simplify(sol_n-true_sphi)==0
# 数値：臨界角より小さい入射では主張の方が sin(phi) を大きく見積もる
V1,V2,W=340.,330.,20.
for sv in np.linspace(0.05,V1/(V2+W)-1e-3,50):
    t=V2*sv/(V1-W*sv); c=(V2+W)*sv/V1
    assert c>t
# 幾何から直接：波面を dt 後に法線方向へ v2 dt 進め，さらに (w dt,0) 平行移動したときの法線方向の変位
ph=np.radians(37.0); n=np.array([np.sin(ph),np.cos(ph)])
disp=V2*n+np.array([W,0.0])
print('normal speed',disp@n,' v2+w sin(phi)=',V2+W*np.sin(ph))
assert abs(disp@n-(V2+W*np.sin(ph)))<1e-9
print('PHY-SND-006 (4) OK')
