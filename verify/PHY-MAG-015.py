# K20092-1 検証：ギャップ幅を考慮したサイクロトロン
import sympy as sp
import numpy as np

q, m, B, V, d, R0, k = sp.symbols('q m B V d R_0 k', positive=True)

# (1) k回目の通過後の速さ・半径
vk = sp.sqrt(2*k*q*V/m)
rk = m*vk/(q*B)
# (2) k回目の通過時間（等加速度 a=qV/(md)）
a = q*V/(m*d)
tau_k = (vk - vk.subs(k, k-1))/a
# 和が d*sqrt(2km/(qV)) になること（望遠鏡和）を k=1..6 で確認
for K in range(1, 7):
    s = sum(tau_k.subs(k, j) for j in range(1, K+1))
    target = d*sp.sqrt(2*K*m/(q*V))
    assert sp.simplify(s - target) == 0
print('gap total ok')

# (3) N, 全時間
N = sp.solve(sp.Eq(rk, R0), k)[0]
print('N =', sp.simplify(N))
t_dee = N*sp.pi*m/(q*B)
t_gap = d*sp.sqrt(2*N*m/(q*V))
print('t_dee =', sp.simplify(t_dee), ' t_gap =', sp.simplify(t_gap))
print('ratio =', sp.simplify(t_gap/t_dee))
K_end = sp.simplify(sp.Rational(1, 2)*m*vk.subs(k, N)**2)
print('K_end =', K_end)

# (4) 固定周期 τ=πm/(qB) で運転：通過完了の累積遅れ <= τ
kmax = sp.solve(sp.Eq(d*sp.sqrt(2*k*m/(q*V)), sp.pi*m/(q*B)), k)[0]
print('kmax =', sp.simplify(kmax))
Vmin = sp.solve(sp.Eq(N, kmax), V)
print('Vmin =', [sp.simplify(x) for x in Vmin])
# t_gap <= tau と同値であることの確認
print('check t_gap/tau at Vmin:', [sp.simplify(t_gap.subs(V, x)/(sp.pi*m/(q*B))) for x in Vmin])

# ---- 数値シミュレーション（同期切り替え）----
qn, mn, Bn, Vn, dn = 1.0, 1.0, 1.0, 0.02, 0.05
Nn = 50
R0n = mn*np.sqrt(2*Nn*qn*Vn/mn)/(qn*Bn)
t = 0.0; v = 0.0; gap_time = 0.0; dee_time = 0.0
an = qn*Vn/(mn*dn)
for j in range(1, Nn+1):
    # ギャップを等加速度で横切る（数値的に細かく積分）
    x = 0.0; dt = 1e-5; vv = v; tt = 0.0
    while x < dn:
        vv += an*dt; x += vv*dt; tt += dt
    gap_time += tt; v = vv
    dee_time += np.pi*mn/(qn*Bn)  # 半円（速さによらない）
print('sim: gap', gap_time, ' formula', dn*np.sqrt(2*Nn*mn/(qn*Vn)))
print('sim: dee', dee_time, ' formula', np.pi*Bn*R0n**2/(2*Vn))
print('sim ratio', gap_time/dee_time, ' formula 2d/(pi R0)', 2*dn/(np.pi*R0n))
# 固定周期運転の数値チェック：累積遅れが τ を超える最初の回
tau = np.pi*mn/(qn*Bn)
lag = 0.0; j = 0; v = 0.0
while True:
    j += 1
    v2 = np.sqrt(2*j*qn*Vn/mn)
    lag += (v2 - v)/an; v = v2
    if lag > tau:
        break
print('first failing crossing', j, ' kmax formula', float(np.pi**2*mn*Vn/(2*qn*Bn**2*dn**2)))
