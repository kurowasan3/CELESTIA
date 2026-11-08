# K20102-1 検証：振動電場中で放された荷電粒子のドリフトと放出時刻
import numpy as np
import sympy as sp
from scipy.integrate import solve_ivp

q, m, E0, T, w, t, t0 = sp.symbols('q m E_0 T omega t t_0', positive=True)

# ---- 方形波 ----
def Esq(tt, E0n, Tn):
    return E0n if (tt % Tn) < Tn/2 else -E0n

qn, mn, E0n, Tn = 1.0, 1.0, 1.0, 1.0
def drift_square(t0n):
    sol = solve_ivp(lambda tt, y: [y[1], qn*Esq(tt, E0n, Tn)/mn], [t0n, t0n+20*Tn], [0, 0],
                    max_step=Tn/2000, rtol=1e-9, atol=1e-12)
    return sol.y[0, -1]/(20*Tn), sol
for t0n in [0.0, 0.1, 0.25, 0.4, 0.6, 0.75, 0.9]:
    vd, _ = drift_square(t0n)
    if t0n < 0.5:
        th = qn*E0n/mn*(Tn/4 - t0n)
    else:
        th = qn*E0n/mn*(t0n - 3*Tn/4)
    print(f'square t0={t0n:.2f} drift sim={vd:+.5f} formula={th:+.5f}')
# t0=T/4 の往復範囲
_, sol = drift_square(0.25)
print('t0=T/4: y range', sol.y[0].min(), sol.y[0].max(), ' formula max', qn*E0n*Tn**2/(16*mn))

# ---- 正弦波 E = E0 cos(wt) ----
v = sp.integrate(q*E0/m*sp.cos(w*t), (t, t0, t))
print('v(t) =', sp.simplify(v))
vd = sp.simplify(v - sp.integrate(v, (t, 0, 2*sp.pi/w))*w/(2*sp.pi) + sp.integrate(v, (t, 0, 2*sp.pi/w))*w/(2*sp.pi))
# 平均速度（1周期平均）
vbar = sp.simplify(sp.integrate(v, (t, t0, t0+2*sp.pi/w))*w/(2*sp.pi))
print('drift =', vbar)
# 時間平均運動エネルギー
Kbar = sp.simplify(sp.integrate(m*v**2/2, (t, t0, t0+2*sp.pi/w))*w/(2*sp.pi))
Up = q**2*E0**2/(4*m*w**2)
print('Kbar/Up =', sp.simplify(Kbar/Up))
# ランダムな放出位相での平均ドリフトエネルギー
Kd = m*vbar**2/2
print('<Kd>/Up =', sp.simplify(sp.integrate(Kd, (t0, 0, 2*sp.pi/w))*w/(2*sp.pi)/Up))
# t0=0 での変位
y = sp.integrate(v.subs(t0, 0), (t, 0, t))
print('y(t) for t0=0:', sp.simplify(y))

# 数値: 正弦波
wn = 2*np.pi
for t0n in [0.0, 0.1, 0.25]:
    sol = solve_ivp(lambda tt, yy: [yy[1], qn*E0n*np.cos(wn*tt)/mn], [t0n, t0n+20], [0, 0], max_step=1e-3, rtol=1e-10, atol=1e-12)
    vdn = sol.y[0, -1]/20
    print(f'sine t0={t0n} drift sim {vdn:+.5f} formula {-qn*E0n/(mn*wn)*np.sin(wn*t0n):+.5f}  ymin {sol.y[0].min():+.4f}')
print('2qE0/(m w^2)=', 2*qn*E0n/(mn*wn**2))
