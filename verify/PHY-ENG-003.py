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
print('constant-force strategy W/(mgH) at that T:', sp.simplify((m*g + 2*m*H/T**2)*H/(m*g*H)).subs(T, Tq))
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
