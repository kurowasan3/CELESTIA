# PHY-MOM-041 : 摩擦の異なる2物体をつなぐひもが反発係数 e で張る運動（相対運動＝床ではずむ球）
import sympy as sp
import numpy as np

e, l, k = sp.symbols('e ell kappa', positive=True)   # kappa = (muB-muA) g cos(theta)
T1 = sp.sqrt(l/k); w1 = k*T1
assert sp.simplify(sp.Rational(1, 2)*k*T1**2 - l/2) == 0
# n回目と n+1 回目の間隔 2 e^n w1 / kappa
n = sp.symbols('n', integer=True, nonnegative=True)
Tinf = T1 + sp.summation(2*e**n*w1/k, (n, 1, sp.oo))
Tinf = sp.simplify(Tinf)
Tinf = Tinf.args[0][0] if isinstance(Tinf, sp.Piecewise) else Tinf
print('T_inf =', sp.simplify(Tinf), ' expected', sp.simplify(T1*(1+e)/(1-e)))
assert sp.simplify(Tinf - T1*(1+e)/(1-e)) == 0
# 損失：換算質量 m/2
m = sp.symbols('m', positive=True)
q = sp.symbols('q', positive=True)  # q = e^2 < 1
geo = sp.summation(q**n, (n, 0, sp.oo))
geo = geo.args[0][0] if isinstance(geo, sp.Piecewise) else geo
loss = sp.Rational(1, 2)*(m/2)*w1**2*(1-q)*geo
loss = sp.simplify(loss)
loss = sp.simplify(loss.args[0][0]) if isinstance(loss, sp.Piecewise) else loss
print('loss =', loss)
assert sp.simplify(loss - m*k*l/4) == 0
# 最小間隔
dmin = l - (e*w1)**2/(2*k)
print('min sep', sp.simplify(dmin))

# ---- 数値シミュレーション（イベント駆動，各区間は等加速度で厳密） ----
def simulate(g=9.8, th=np.radians(30), muA=0.1, muB=0.4, L=1.0, v0=2.0, ee=0.6, mm=1.0, tmax=10):
    aA = g*(np.sin(th)-muA*np.cos(th)); aB = g*(np.sin(th)-muB*np.cos(th))
    kap = aA - aB
    t = 0.0; xA = L/2; xB = 0.0; vA = v0; vB = v0
    taut_times = []; Eloss = 0.0; minsep = 1e9; minvA = 1e9; minvB = 1e9
    while t < tmax and len(taut_times) < 200:
        # separation s = xA-xB, ds = vA-vB, dd = kap
        s0 = xA-xB; w = vA-vB
        # time until s = L : s0 + w t + kap t^2/2 = L
        disc = w*w + 2*kap*(L-s0)
        dt = (-w + np.sqrt(disc))/kap
        # sample min separation & velocities in interval
        ts = np.linspace(0, dt, 200)
        minsep = min(minsep, np.min(s0 + w*ts + kap*ts**2/2))
        minvA = min(minvA, np.min(vA + aA*ts)); minvB = min(minvB, np.min(vB + aB*ts))
        xA += vA*dt + aA*dt**2/2; xB += vB*dt + aB*dt**2/2
        vA += aA*dt; vB += aB*dt; t += dt
        taut_times.append(t)
        V = (vA+vB)/2; wr = vA-vB
        KE0 = 0.5*mm*(vA**2+vB**2)
        vA = V - ee*wr/2; vB = V + ee*wr/2
        Eloss += KE0 - 0.5*mm*(vA**2+vB**2)
        if wr < 1e-12:
            break
    return dict(kap=kap, taut=taut_times, Eloss=Eloss, minsep=minsep, minvA=minvA, minvB=minvB,
                g=g, th=th, muA=muA, muB=muB, L=L, ee=ee, mm=mm)
r = simulate()
kap = r['kap']; L = 1.0; ee = 0.6
T1n = np.sqrt(L/kap)
print('first taut', r['taut'][0], 'theory', T1n)
print('limit', r['taut'][-1], 'theory', T1n*(1+ee)/(1-ee))
print('Eloss', r['Eloss'], 'theory', 1.0*kap*L/4)
print('min sep', r['minsep'], ' min vA', r['minvA'], ' min vB', r['minvB'])
assert abs(r['taut'][0]-T1n) < 1e-9
assert abs(r['taut'][-1]-T1n*(1+ee)/(1-ee)) < 1e-6
assert abs(r['Eloss']-kap*L/4) < 1e-8
assert r['minsep'] >= L/2 - 1e-9 and r['minvA'] > 0 and r['minvB'] > 0
# 間隔の比 e
d = np.diff(r['taut'])
print('interval ratios', d[1:5]/d[:4])
# e=1：間隔一定 2T1
r1 = simulate(ee=1.0, tmax=5)
d1 = np.diff(r1['taut'])
print('e=1 intervals', d1[:4], ' 2T1=', 2*T1n, ' loss', r1['Eloss'])
assert np.allclose(d1, 2*T1n) and abs(r1['Eloss']) < 1e-9
# 損失が e によらないことの確認
for ee2 in [0.2, 0.5, 0.9]:
    rr = simulate(ee=ee2, tmax=60)
    assert abs(rr['Eloss']-kap*L/4) < 1e-6, (ee2, rr['Eloss'])
# 一体運動のときの張力  A: m a = m g sin - muA m g cos - T, B: m a = m g sin - muB m g cos + T
g = 9.8; th = np.radians(30)
T = 0.5*(0.4-0.1)*g*np.cos(th)
ac = g*np.sin(th) - 0.5*(0.1+0.4)*g*np.cos(th)
assert abs((g*np.sin(th) - 0.1*g*np.cos(th) - T) - ac) < 1e-12 and abs((g*np.sin(th) - 0.4*g*np.cos(th) + T) - ac) < 1e-12
# v0 条件  v0 > w1/2
print('w1/2 =', np.sqrt(kap*L)/2)
print('ALL OK')
