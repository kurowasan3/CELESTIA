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
# 数値
V1,V2,W,H=340.,330.,20.,300.
sc=V1/(V2+W); thc=np.arcsin(sc)
print('theta_c deg',np.degrees(thc),' skip distance 2H tan',2*H*np.tan(thc))
su=V1*W/(V2**2-W**2); print('upwind vertical theta deg',np.degrees(np.arcsin(su)), su)
# 直接：風上で sin th = su のとき sin ph = W/V2 となるか
print(V2*su/(V1+W*su), W/V2)
# 全反射の直前で光線の向き（ph->90°）: ray = (v2+w, 0) 水平
