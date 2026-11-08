# K20121-1 検証：粗い床上で打ち出された板と上の小物体（停止時刻が μ2 によらない条件）
import sympy as sp
import numpy as np

M, m, g, v0, mu1, mu2, t = sp.symbols('M m g v_0 mu_1 mu_2 t', positive=True)

aB = -(mu1*(M+m) + mu2*m)*g/M   # 板の加速度（初め）
ab = mu2*g                        # 小物体
t1 = sp.solve(sp.Eq(v0 + aB*t, ab*t), t)[0]
v1 = sp.simplify(ab*t1)
print('t1 =', sp.simplify(t1), ' v1 =', v1)
T_mom = M*v0/(mu1*(M+m)*g)
# μ1<=μ2 : 一体で減速 mu1 g → 停止 t1 + v1/(mu1 g)
T_tog = sp.simplify(t1 + v1/(mu1*g))
print('together stop time - momentum formula =', sp.simplify(T_tog - T_mom))
# μ1>μ2 : 再び滑る
aB2 = -(mu1*(M+m) - mu2*m)*g/M
tb = sp.simplify(t1 + v1/(-aB2))
ts = sp.simplify(t1 + v1/(mu2*g))
print('board stops tb =', tb)
print('block stops ts =', ts, '  ts/T_mom =', sp.simplify(ts/T_mom))
print('tb < ts  <=>  (aB2 magnitude - mu2 g) =', sp.factor(sp.simplify(-aB2 - mu2*g)))

# 数値シミュレーション（摩擦の場合分けを逐次判定）
def simulate(Mn, mn, mu1n, mu2n, v0n, gn=9.8, dt=1e-6):
    VB, vb = v0n, 0.0
    t_ = 0.0; tB = None
    while True:
        rel = vb - VB
        # 板と物体の間の摩擦力（物体が受ける）
        if abs(rel) > 1e-9:
            f2 = -np.sign(rel)*mu2n*mn*gn
        else:
            # 一体で動けるか判定
            if VB > 1e-12:
                a_common = -mu1n*gn
                need = mn*a_common
                if abs(need) <= mu2n*mn*gn:
                    f2 = need             # 一体で減速できる（静止摩擦）
                else:
                    f2 = -mu2n*mn*gn      # 物体は前へ滑り出す（摩擦は後ろ向き）
            else:
                f2 = -np.sign(vb)*mu2n*mn*gn if abs(vb) > 1e-12 else 0.0
        # 床の摩擦（板が受ける）
        N1 = (Mn+mn)*gn
        if VB > 1e-12:
            f1 = -mu1n*N1
        else:
            # 板は静止: 物体からの力 -f2 を静止摩擦で支えられるか
            if abs(-f2) <= mu1n*N1:
                f1 = f2   # つり合う
            else:
                f1 = -np.sign(-f2)*mu1n*N1
            if tB is None:
                tB = t_
        aB_ = (f1 - f2)/Mn
        ab_ = f2/mn
        VB_new = VB + aB_*dt; vb_new = vb + ab_*dt
        if VB > 0 and VB_new <= 0: VB_new = 0.0
        if vb > 0 and vb_new <= 0 and abs(VB) < 1e-12: vb_new = 0.0
        # 相対速度の符号が変わったらそろえる
        if (vb - VB)*(vb_new - VB_new) < 0:
            vc = (mn*vb_new + Mn*VB_new)/(mn+Mn); vb_new = VB_new = vc
        VB, vb = VB_new, vb_new
        t_ += dt
        if VB <= 0 and vb <= 0:
            return t_, tB
        if t_ > 10: return None

Mn, mn, v0n, gn = 2.0, 1.0, 3.0, 9.8
for mu1n, mu2n in [(0.2, 0.5), (0.4, 0.2)]:
    res = simulate(Mn, mn, mu1n, mu2n, v0n)
    sub = {M: Mn, m: mn, mu1: mu1n, mu2: mu2n, v0: v0n, g: gn}
    print(f'mu1={mu1n}, mu2={mu2n}: sim all-stop {res[0]:.5f}, board stop {res[1]}, '
          f'T_mom={float(T_mom.subs(sub)):.5f}, ts={float(ts.subs(sub)):.5f}, tb={float(tb.subs(sub)):.5f}')
