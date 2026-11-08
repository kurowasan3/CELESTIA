import numpy as np
from scipy.optimize import minimize_scalar, brentq
c=340.0; cR=350.0; L=10.0; D=5.0; f=500.0; W=1.7
def tau(y):
    r=minimize_scalar(lambda b: np.hypot(L,b)/c+np.hypot(D,y-b)/cR, bounds=(-1e3,1e3), method='bounded', options={'xatol':1e-12})
    return r.fun, r.x
for y in (2.0,8.775,30.0,300.0):
    t,b=tau(y); h=1e-4
    dtdy=(tau(y+h)[0]-tau(y-h)[0])/(2*h)
    fQ=f*(1+W*dtdy)   # Q moves -y: received phase rate
    sth=b/np.hypot(L,b); sph=(y-b)/np.hypot(D,y-b)
    print(f'y={y}: fQ={fQ:.5f}, f(1+W sin th/c)={f*(1+W*sth/c):.5f}, f(1+W sin ph/cR)={f*(1+W*sph/cR):.5f}, snell {sth/c:.6e} {sph/cR:.6e}')
print('limit fmax', f*(1+W/cR), ' critical sin', c/cR)
# (5) 新設問：到達時間 tau(y) によるドップラー効果
import sympy as sp
bs,ys,Ls,Ds,cs,cRs=sp.symbols('b y L D c c_R',positive=True)
tt=sp.sqrt(Ls**2+bs**2)/cs+sp.sqrt(Ds**2+(ys-bs)**2)/cRs
dtdb=sp.diff(tt,bs)
sin_th=bs/sp.sqrt(Ls**2+bs**2); sin_ph=(ys-bs)/sp.sqrt(Ds**2+(ys-bs)**2)
assert sp.simplify(dtdb-(sin_th/cs-sin_ph/cRs))==0
d2=Ls**2/(cs*(Ls**2+bs**2)**sp.Rational(3,2))+Ds**2/(cRs*(Ds**2+(ys-bs)**2)**sp.Rational(3,2))
assert sp.simplify(sp.diff(dtdb,bs)-d2)==0   # 正の量の和 → dt/db は b の増加関数（極小）
print('d2t/db2 = L^2/(c(L^2+b^2)^1.5) + D^2/(c_R(...)^1.5) > 0')
assert sp.simplify(sp.diff(tt,ys)-sin_ph/cRs)==0
# 数値：dtau/dy（最小化した到達時間の数値微分）= sin(theta)/c，f_Q の一致
for yv in (1.0,5.0,20.0,200.0):
    t0,b0=tau(yv); h=1e-4
    d=(tau(yv+h)[0]-tau(yv-h)[0])/(2*h)
    sth=b0/np.hypot(L,b0)
    print(f'y={yv}: dtau/dy={d:.8e}, sin(th)/c={sth/c:.8e}')
    assert abs(d-sth/c)<1e-7
# 直接シミュレーション：Q を動かし，到達する音の発射時刻の変化率から受信振動数を求める
def emit_time(tq):
    yq=50.0-W*tq; return tq-tau(yq)[0]
tq=3.0; h=1e-3
fQ_sim=f*(emit_time(tq+h)-emit_time(tq-h))/(2*h)
b0=tau(50.0-W*tq)[1]; sth=b0/np.hypot(L,b0)
print('f_Q simulated',fQ_sim,' formula',f*(1+W*sth/c))
assert abs(fQ_sim-f*(1+W*sth/c))<1e-4
print('PHY-SND-010 (5) OK')
