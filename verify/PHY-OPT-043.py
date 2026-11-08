# PHY-OPT-043 : 位相がゆらぐ二重スリット干渉の強度相関（一様ゆらぎと二値ゆらぎ）
import sympy as sp
import numpy as np

A, Bv, ph = sp.symbols('A B phi', real=True)
al = sp.symbols('alpha', positive=True)
I1 = al/2*(1 + sp.cos(A - ph)); I2 = al/2*(1 + sp.cos(Bv - ph))
# 一様分布の平均
avg_u = lambda f: sp.simplify(sp.integrate(f, (ph, 0, 2*sp.pi))/(2*sp.pi))
Ibar = avg_u(I1)
Mbar = sp.simplify(sp.expand_trig(avg_u(sp.expand(I1*I2))))
print('uniform: Ibar', Ibar, ' Mbar', sp.simplify(Mbar))
assert sp.simplify(Ibar - al/2) == 0
assert sp.simplify(Mbar - al**2/4*(1 + sp.cos(A - Bv)/2)) == 0
# 二値 φ = 0, π
avg_b = lambda f: sp.simplify((f.subs(ph, 0) + f.subs(ph, sp.pi))/2)
Ib = avg_b(I1); Mb = avg_b(sp.expand(I1*I2))
print('binary: Ibar', Ib, ' Mbar', sp.simplify(Mb))
assert sp.simplify(Ib - al/2) == 0
assert sp.simplify(Mb - al**2/4*(1 + sp.cos(A)*sp.cos(Bv))) == 0

# 幾何：B が +x 側，B の光が φ 遅れる → 位相差 δ = (2π/λ)(L_A - L_B) - φ,  L_A - L_B ≈ a x / L
lam, L, a, x = sp.symbols('lambda L a x', positive=True)
LA = sp.sqrt(L**2 + (x + a/2)**2); LB = sp.sqrt(L**2 + (x - a/2)**2)
eps = sp.symbols('eps', positive=True)
dexpr = (LA - LB).subs({x: eps*x, a: eps*a}, simultaneous=True)
diff = sp.series(dexpr, eps, 0, 3).removeO().subs(eps, 1)
print('L_A - L_B ≈', sp.simplify(diff))
assert sp.simplify(diff - a*x/L) == 0
# (4) 数値
P = 0.40e-3; lamv = 600e-9; Lv = 2.0
print('a =', lamv*Lv/P)
# モンテカルロで (3)(5) を確認
rng = np.random.default_rng(1)
Pn = 1.0; x1 = 0.13; x2s = np.linspace(-1, 1, 9)
phs = rng.uniform(0, 2*np.pi, 400000)
for x2 in x2s:
    r1 = 0.5*(1+np.cos(2*np.pi*x1/Pn - phs)); r2 = 0.5*(1+np.cos(2*np.pi*x2/Pn - phs))
    mc = np.mean(r1*r2); th = 0.25*(1+0.5*np.cos(2*np.pi*(x2-x1)/Pn))
    assert abs(mc-th) < 3e-3, (x2, mc, th)
phb = rng.choice([0, np.pi], 400000)
for x2 in x2s:
    r1 = 0.5*(1+np.cos(2*np.pi*x1/Pn - phb)); r2 = 0.5*(1+np.cos(2*np.pi*x2/Pn - phb))
    mc = np.mean(r1*r2); th = 0.25*(1+np.cos(2*np.pi*x1/Pn)*np.cos(2*np.pi*x2/Pn))
    assert abs(mc-th) < 3e-3, (x2, mc, th)
# 二値で x1 = P/4 なら x2 によらない
x1 = 0.25
vals = [0.25*(1+np.cos(2*np.pi*x1)*np.cos(2*np.pi*x2)) for x2 in x2s]
assert max(vals)-min(vals) < 1e-12
# (6)(a) same position: g-1 = var(I)/mean(I)^2 >= 0 for arbitrary phi distributions
for trial in range(200):
    k = rng.integers(1, 6); vals_phi = rng.uniform(0, 2*np.pi, k); w = rng.uniform(0, 1, k); w /= w.sum()
    x = rng.uniform(0, 1)
    I = 0.5*(1+np.cos(2*np.pi*x - vals_phi))
    m = np.sum(w*I); M = np.sum(w*I*I); var = np.sum(w*(I-m)**2)
    assert abs((M/m**2 - 1) - var/m**2) < 1e-12 and M/m**2 >= 1-1e-12
print('(6)(a) g(x,x) >= 1 for random discrete phi distributions; equality iff var=0')
# binary at x=P/4: var = 0 -> g = 1
I = 0.5*(1+np.cos(2*np.pi*0.25 - np.array([0, np.pi]))); assert np.ptp(I) < 1e-12
# (6)(b) binary case: g = 0 exactly at x1 = jP/2, x2 = x1 + (k+1/2)P
gb = lambda x1, x2: 1+np.cos(2*np.pi*x1)*np.cos(2*np.pi*x2)
X1, X2 = np.meshgrid(np.linspace(-1, 1, 801), np.linspace(-1, 1, 801))
G = gb(X1, X2); zs = np.argwhere(G < 1e-12)
for (r, c) in zs:
    x1v, x2v = X1[r, c], X2[r, c]
    assert abs(2*x1v - round(2*x1v)) < 1e-9 and abs((x2v-x1v-0.5) - round(x2v-x1v-0.5)) < 1e-9
print('(6)(b) zeros of binary g found at', len(zs), 'grid points, all of the predicted form')
for (x1v, x2v) in [(0, 0.5), (0.5, 0), (0.5, 1.0), (1.5, -1.0)]:
    r1 = 0.5*(1+np.cos(2*np.pi*x1v - np.array([0, np.pi]))); r2 = 0.5*(1+np.cos(2*np.pi*x2v - np.array([0, np.pi])))
    assert np.allclose(r1*r2, 0)
# uniform case never reaches 0: min g = 1/2
print('uniform g min =', min(1+0.5*np.cos(2*np.pi*d) for d in np.linspace(0, 1, 1001)))
print('ALL OK')
