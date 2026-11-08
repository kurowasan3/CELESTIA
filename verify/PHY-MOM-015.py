# K20111-1 検証：砂を集めながら一定の力で押される板（可変質量と摩擦）
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp

M, sig, mu, g, F, u, s, m, X = sp.symbols('M sigma mu g F u s m X', positive=True)

# (1) 等速 u: 力 = mu g (M+sig s) + sig u^2, 取り込みで失われる仕事率 sig u^3/2
Fu = mu*g*(M+sig*s) + sig*u**2
P_in = Fu*u
P_fric = mu*g*(M+sig*s)*u
dK = sig*u*u**2/2          # 単位時間に取り込まれる砂の運動エネルギー
print('capture loss rate =', sp.simplify(P_in - P_fric - dK))

# (2) 一定の力 F: (m v)^2 = [F(m^2-M^2) - 2/3 mu g (m^3-M^3)]/sig
P2 = (F*(m**2-M**2) - sp.Rational(2, 3)*mu*g*(m**3-M**3))/sig
# 微分方程式 d(P^2/2)/dm = (F - mu g m) m / sig を満たすか
print('ODE check:', sp.simplify(sp.diff(P2/2, m) - (F - mu*g*m)*m/sig), ' P2(M)=', sp.simplify(P2.subs(m, M)))

# (3) F = 2 mu M g : 停止位置
P2X = sp.simplify(P2.subs({F: 2*mu*M*g, m: X*M}))
sols = sp.solve(sp.Eq(P2X, 0), X)
print('stop X:', sols)
# (4) 最大速さ: v^2 = P2/m^2
v2 = sp.simplify(P2X/(X*M)**2)
Xm = sp.solve(sp.diff(v2, X), X)
print('Xmax:', Xm, ' v_max^2 =', sp.simplify(v2.subs(X, Xm[0])))
# (5) エネルギー: 仕事, 床の摩擦熱, 取り込みによる損失
Xs = 1 + sp.sqrt(3)
W = 2*mu*M*g*(Xs-1)*M/sig
Q_fric = mu*g*((Xs*M)**2 - M**2)/(2*sig)
loss = sp.simplify(W - Q_fric)
print('W =', sp.simplify(W), ' Q_fric =', sp.simplify(Q_fric), ' loss =', loss, ' loss/W =', sp.simplify(loss/W), float(loss/W))

# ---- 数値積分（時間積分）----
Mn, sn, mun, gn = 1.0, 2.0, 0.3, 9.8
Fn = 2*mun*Mn*gn
def rhs(t, y):
    x, P = y           # P = (M + sigma x) v
    mm = Mn + sn*x
    v = P/mm
    return [v, Fn - mun*gn*mm]
def stop(t, y): return y[1]
stop.terminal = True; stop.direction = -1
sol = solve_ivp(rhs, [0, 50], [0, 0], events=stop, max_step=1e-3, rtol=1e-10, atol=1e-12)
xs = sol.y[0, -1]
print('sim stop s =', xs, ' formula sqrt3 M/sigma =', np.sqrt(3)*Mn/sn)
v = sol.y[1]/(Mn + sn*sol.y[0])
i = np.argmax(v)
print('sim vmax =', v[i], ' at m/M =', 1 + sn*sol.y[0, i]/Mn, ' formula', np.sqrt((2-4**(1/3))*mun*gn*Mn/sn), 4**(1/3))
# 摩擦熱の数値積分
Qf = np.trapezoid(mun*gn*(Mn + sn*sol.y[0]), sol.y[0])
Wn = Fn*xs
print('sim loss/W =', (Wn - Qf)/Wn, ' formula', 0.5 - np.sqrt(3)/4)
# 停止後: 押す力 F と最大摩擦力の比較
print('F =', Fn, ' max static friction at stop =', mun*gn*(Mn + sn*xs))
