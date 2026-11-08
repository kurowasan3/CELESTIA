# K20093-1 検証：浅い水のくさび形波面と深い水のケルビン角
import sympy as sp
import numpy as np

V, c, g, k, dk, t, phi = sp.symbols('V c g k Delta_k t phi', positive=True)

# (1) 浅い水: sin(theta)=c/V （素元波の包絡線）  数値で確認
Vn, cn = 2.0, 1.0
ang = []
for tau in np.linspace(0.01, 1, 200):          # tau 前に P=(−V tau,0) で出た素元波, 船は原点
    # 船(原点)からこの円への接線の角 = asin(c tau / (V tau))
    ang.append(np.arcsin(cn*tau/(Vn*tau)))
print('shallow: half angle', np.degrees(max(ang)), ' asin(c/V)=', np.degrees(np.arcsin(cn/Vn)))

# (2) 群速度: omega = sqrt(g k)
w = sp.sqrt(g*k)
vg = sp.limit((w.subs(k, k+dk) - w.subs(k, k-dk))/(2*dk), dk, 0)
cph = w/k
print('vg/c =', sp.simplify(vg/cph))

# 数値: 2波の重ね合わせの包絡線の移動速度
gn = 9.8; k0 = 1.0; d = 0.02
w1 = np.sqrt(gn*(k0+d)); w2 = np.sqrt(gn*(k0-d))
print('numeric envelope speed', (w1-w2)/(2*d), ' c/2 =', np.sqrt(gn/k0)/2)

# (3) 定常条件 c = V cos(phi) → 波長
lam = 2*sp.pi*(V*sp.cos(phi))**2/g
print('lambda(phi) =', lam)
# check: c(lam)=sqrt(g lam/2pi)
print('check c:', sp.simplify(sp.sqrt(g*lam/(2*sp.pi)) - V*sp.cos(phi)))

# (4) エネルギーの位置 Q = P + (V t cos(phi)/2)(cos phi, sin phi), P=原点 → 円
X = V*t*sp.cos(phi)**2/2; Y = V*t*sp.cos(phi)*sp.sin(phi)/2
print('circle check:', sp.simplify((X - V*t/4)**2 + Y**2 - (V*t/4)**2))

# (5) 船 S=(V t,0) から見た最大角 → asin(1/3)
Vn = 1.0
best = 0; bestphi = None
for tt in [0.5, 1.0, 3.0]:
    ph = np.linspace(-np.pi/2+1e-6, np.pi/2-1e-6, 200001)
    Xn = Vn*tt*np.cos(ph)**2/2; Yn = Vn*tt*np.cos(ph)*np.sin(ph)/2
    a = np.arctan2(np.abs(Yn), Vn*tt - Xn)
    i = np.argmax(a)
    if a[i] > best:
        best = a[i]; bestphi = ph[i]
print('wake half angle (deg)', np.degrees(best), ' asin(1/3)=', np.degrees(np.arcsin(1/3)))
print('phi at edge (deg)', np.degrees(abs(bestphi)), ' atan(1/sqrt2)=', np.degrees(np.arctan(1/np.sqrt(2))))
print('cos^2 phi_edge', np.cos(bestphi)**2, ' lambda_edge factor (×V^2/g)', 2*np.pi*np.cos(bestphi)**2, 4*np.pi/3)
# 境界上の点 T の S からの距離 / (V t) = 1/sqrt2
