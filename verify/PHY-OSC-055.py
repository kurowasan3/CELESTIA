# PHY-OSC-055 : 回転する水平管から飛び出す小球（回転系と床の系の乗りかえ）
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp

w, L, r0, tau, m, x = sp.symbols('omega L r_0 tau m x', positive=True)
# (1) 回転系のエネルギー保存（遠心力の位置エネルギー -1/2 m w^2 x^2）
u = w*sp.sqrt(L**2 - r0**2)
# 独立に：r'' = w^2 r の解 r = r0 cosh(w t)
t = sp.symbols('t', positive=True)
rt = r0*sp.cosh(w*t)
assert sp.simplify(sp.diff(rt, t, 2) - w**2*rt) == 0
assert sp.simplify(sp.diff(rt, t)**2 - w**2*(rt**2 - r0**2)) == 0
# (2) 側壁の力：床の系で角運動量 m w x^2 の変化率 = N x  -> N = 2 m w xdot
N = sp.diff(m*w*rt**2, t)/rt
assert sp.simplify(N - 2*m*w*sp.diff(rt, t)) == 0
# (3) 落下点
X = L + u*tau; Y = w*L*tau
D = sp.sqrt(X**2 + Y**2)
print('D =', sp.simplify(D))
# (4) 45°条件  w L tau = L + u tau
r0sol = sp.solve(sp.Eq(w*L*tau, L + u*tau), r0)
print('r0 for 45deg:', [sp.simplify(z) for z in r0sol])
r0_formula = L*sp.sqrt(1 - (1 - 1/(w*tau))**2)
for z in r0sol:
    print(sp.simplify(z**2 - r0_formula**2))
# (5) 数値： r0 = sqrt3/2 L, w tau = 1
num = {L: 1, w: 1, tau: 1, r0: sp.sqrt(3)/2}
psi = 1 - sp.atan(sp.Rational(2, 3))
print('psi rad', float(psi), 'deg', float(psi*180/sp.pi), 'D/L', float(D.subs(num)))

# solve_ivp による独立確認（床の系で，管の中は拘束運動 → 出口以降は放物運動）
def run(omega, Lv, r0v, h, g=9.8):
    def f(tt, y):  # 回転系の動径運動
        return [y[1], omega**2*y[0]]
    ev = lambda tt, y: y[0]-Lv
    ev.terminal = True
    s = solve_ivp(f, [0, 50], [r0v, 0], events=ev, rtol=1e-11, atol=1e-13)
    te = s.t_events[0][0]; ur = s.y_events[0][0][1]
    th_e = omega*te   # 出口での管の角度
    # 床の系の速度
    er = np.array([np.cos(th_e), np.sin(th_e)]); et = np.array([-np.sin(th_e), np.cos(th_e)])
    pos = Lv*er; vel = ur*er + omega*Lv*et
    tf = np.sqrt(2*h/g)
    land = pos + vel*tf
    tube_angle = th_e + omega*tf
    land_angle = np.arctan2(land[1], land[0])
    psi_v = (tube_angle - land_angle)
    return ur, np.linalg.norm(land), psi_v, land_angle - th_e
g = 9.8; omega = 2.0; Lv = 0.5
h = g/(2*omega**2)   # omega*tau = 1
ur, Dn, psin, phin = run(omega, Lv, np.sqrt(3)/2*Lv, h)
print('num u', ur, 'theory', omega*Lv/2, 'D/L', Dn/Lv, 'psi', psin)
assert abs(ur-omega*Lv/2) < 1e-7 and abs(Dn/Lv-np.sqrt(13)/2) < 1e-6 and abs(psin-float(psi)) < 1e-6
# 45° 条件の数値確認： omega*tau = 2
h2 = g*4/(2*omega**2)
r0v = Lv*np.sqrt(1-(1-1/2)**2)
_, _, _, phi2 = run(omega, Lv, r0v, h2)
print('phi for 45 cond', np.degrees(phi2))
assert abs(np.degrees(phi2)-45) < 1e-5
# psi>0 をランダムに確認
rng = np.random.default_rng(0)
for _ in range(2000):
    om = rng.uniform(0.2, 5); Lr = rng.uniform(0.1, 2); r0r = rng.uniform(0.01, 0.999)*Lr; hh = rng.uniform(0.05, 3)
    assert run(om, Lr, r0r, hh)[2] > 0
print('ALL OK')
