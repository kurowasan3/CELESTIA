# T20081-1 minimum work to lift to height H passing at time T, 0<=F<=F0
import sympy as sp
import numpy as np
from scipy.optimize import minimize

m, g, H, F0, T, tau, F1 = sp.symbols('m g H F_0 T tau F_1', positive=True)
# (1) T_min
Tmin = sp.sqrt(2*m*H/(F0 - m*g))
print('Tmin', Tmin)
# (3) T*: bang F0 for tau then free flight, apex at H exactly at time T
a0 = F0/m - g
eqs = [sp.Eq(a0*tau, g*(T - tau)), sp.Eq(a0*tau**2/2 + (a0*tau)**2/(2*g), H)]
sol = sp.solve(eqs, [tau, T], dict=True)
print('T* candidates', [sp.simplify(s_[T]) for s_ in sol])
Tstar = sp.sqrt(2*H*F0/(g*(F0 - m*g)))
print('T* claim check', [sp.simplify(s_[T] - Tstar) for s_ in sol])
# (4) F0=2mg, Tmin<T<T*: bang tau then coast, reach H at T
u = sp.sqrt(T**2/2 - H/g)
tau4 = T - u
xT = g*tau4**2/2 + g*tau4*(T - tau4) - g*(T - tau4)**2/2
print('x(T)-H', sp.simplify(xT - H))
W4 = 2*m*g*(g*tau4**2/2)
vT = g*tau4 - g*(T - tau4)
print('W - (mgH + m vT^2/2)', sp.simplify(sp.expand(W4 - m*g*H - m*vT**2/2)))
Tq = sp.sqrt(3*H/g)
print('W at T=sqrt(3H/g):', sp.nsimplify(sp.simplify(W4.subs(T, Tq)/(m*g*H))), float(W4.subs(T, Tq)/(m*g*H)))
print('v_T at that T /sqrt(gH):', sp.simplify(vT.subs(T, Tq)/sp.sqrt(g*H)))

# Numerical optimisation over piecewise-constant force profiles (m=g=H=1, F0=2)
def sim(Fs, Ttot):
    N = len(Fs); dt = Ttot/N; x = 0.0; v = 0.0; W = 0.0
    for F in Fs:
        a = F - 1.0
        if x <= 0 and v <= 0 and a < 0:
            a = 0.0  # resting on floor
        x_new = x + v*dt + a*dt*dt/2
        W += F*(x_new - x)
        v += a*dt; x = x_new
    return x, v, W

for Ttot in [np.sqrt(3.0), 1.6, 2.2]:
    N = 60
    cons = [{'type': 'eq', 'fun': lambda Fs: sim(Fs, Ttot)[0] - 1.0}]
    best = None
    for trial in range(6):
        x0 = np.clip(np.random.rand(N)*2, 0, 2)
        r = minimize(lambda Fs: sim(Fs, Ttot)[2], x0, bounds=[(0, 2)]*N, constraints=cons, method='SLSQP', options={'maxiter': 500})
        if r.success and (best is None or r.fun < best.fun):
            best = r
    Ts = 2.0
    if Ttot >= Ts:
        theory = 1.0
    else:
        uu = np.sqrt(Ttot**2/2 - 1); theory = (Ttot - uu)**2*1.0
    print('T=%.4f  numeric min W=%.4f  theory=%.4f' % (Ttot, best.fun, theory))

# ---- (5) general F0: bang F0 for tau, then coast ----
eps = m*(g*T**2 + 2*H)/F0
tauG = T - sp.sqrt(T**2 - eps)
a0G = F0/m - g
xTG = a0G*tauG**2/2 + a0G*tauG*(T - tauG) - g*(T - tauG)**2/2
print('(5) x(T)-H:', sp.simplify(sp.expand(xTG - H)))
WG = F0*a0G*tauG**2/2
vTG = a0G*tauG - g*(T - tauG)
print('(5) W - (mgH + m vT^2/2):', sp.simplify(sp.expand(WG - m*g*H - m*vTG**2/2)))
print('(5) F0=2mg reduces to (4):', sp.simplify(WG.subs(F0, 2*m*g) - W4))
Wlim = sp.limit(WG, F0, sp.oo)
print('(5) F0->oo limit:', sp.factor(Wlim), ' claim diff:', sp.simplify(Wlim - m*(H/T + g*T/2)**2/2))
# monotonic decrease in F0 toward the limit (numeric sample m=g=H=1, T=1.2)
vals = [float(WG.subs({m:1, g:1, H:1, T:1.2, F0:f})) for f in [3, 5, 10, 100, 1e4]]
print('(5) W(F0) samples:', vals, ' limit', float(Wlim.subs({m:1, g:1, H:1, T:1.2})))
assert all(vals[i] > vals[i+1] for i in range(len(vals)-1))
# numeric optimisation with F0=3 (m=g=H=1), T=1.2 (Tmin=1, T*=sqrt(3))
def sim3(Fs, Ttot):
    N = len(Fs); dt = Ttot/N; x = 0.0; v = 0.0; W = 0.0
    for F in Fs:
        a = F - 1.0
        if x <= 0 and v <= 0 and a < 0:
            a = 0.0
        x_new = x + v*dt + a*dt*dt/2
        W += F*(x_new - x); v += a*dt; x = x_new
    return x, v, W
Ttot = 1.2; N = 60; best = None
np.random.seed(1)
for trial in range(6):
    x0 = np.random.rand(N)*3
    r = minimize(lambda Fs: sim3(Fs, Ttot)[2], x0, bounds=[(0, 3)]*N,
                 constraints=[{'type': 'eq', 'fun': lambda Fs: sim3(Fs, Ttot)[0] - 1.0}], method='SLSQP', options={'maxiter': 500})
    if r.success and (best is None or r.fun < best.fun):
        best = r
th = float(WG.subs({m:1, g:1, H:1, T:Ttot, F0:3}))
print('(5) F0=3,T=1.2 numeric min W=%.4f theory=%.4f' % (best.fun, th))
assert abs(best.fun - th) < 0.02
