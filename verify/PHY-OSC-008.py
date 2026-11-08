# M33-1 verification: ball on spring-supported plate, elastic, time-reversal "return" condition
import sympy as sp
import numpy as np
from scipy.optimize import brentq

m, M, k, g, h = sp.symbols('m M k g h', positive=True)
v0 = sp.sqrt(2*g*h)
v = (m-M)/(m+M)*v0          # ball after (down positive)
V = 2*m/(m+M)*v0            # plate after
# check momentum and e=1
print('momentum:', sp.simplify(m*v + M*V - m*v0), ' e:', sp.simplify((V - v)/v0))
A = V*sp.sqrt(M/k)
print('(2) amplitude =', sp.simplify(A))
cond = sp.Eq(2*(-v)/g, sp.pi*sp.sqrt(M/k))
hsol = sp.solve(cond, h)
print('(3) h =', [sp.factor(s) for s in hsol])
# second collision: ball +|v| (down), plate -V (up)
vb, vp = sp.symbols('vb vp')
s2 = sp.solve([sp.Eq(m*vb + M*vp, m*(-v) + M*(-V)), sp.Eq(vp - vb, -((-V) - (-v)))], [vb, vp], dict=True)[0]
print('(4) after 2nd: ball', sp.simplify(s2[vb]/v0), 'v0 ; plate', sp.simplify(s2[vp]))

# ---------- numerical event simulation (down positive, origin = plate equilibrium = collision point) ----------
def simulate(mv, Mv, kv, gv, hv, tmax_factor=3.0):
    w = np.sqrt(kv/Mv)
    # state: ball position yb (down +), vb ; plate xp, vp. Ball starts at -hv at rest.
    t = 0.0; yb, vb = -hv, 0.0; xp, vp = 0.0, 0.0
    log = []
    T_end = tmax_factor*(2*np.sqrt(2*hv/gv) + np.pi/w)
    while t < T_end and len(log) < 6:
        def gap(tau):   # plate position - ball position (>0 : ball above plate)
            yb_t = yb + vb*tau + 0.5*gv*tau**2
            xp_t = xp*np.cos(w*tau) + vp/w*np.sin(w*tau)
            return xp_t - yb_t
        # find first root after small epsilon
        taus = np.linspace(1e-9, T_end, 400001)
        gv_ = np.array([gap(s) for s in taus[:1]])
        vals = gap(taus)
        idx = np.where(np.sign(vals[:-1]) != np.sign(vals[1:]))[0]
        # also tangential/near-zero minima
        if len(idx) == 0:
            break
        i = idx[0]
        tau = brentq(gap, taus[i], taus[i+1], xtol=1e-14)
        yb_c = yb + vb*tau + 0.5*gv*tau**2; vb_c = vb + gv*tau
        xp_c = xp*np.cos(w*tau) + vp/w*np.sin(w*tau); vp_c = -xp*w*np.sin(w*tau) + vp*np.cos(w*tau)
        # elastic collision
        nb = ((mv-Mv)*vb_c + 2*Mv*vp_c)/(mv+Mv)
        np_ = ((Mv-mv)*vp_c + 2*mv*vb_c)/(mv+Mv)
        t += tau
        log.append((t, xp_c, vb_c, vp_c, nb, np_))
        yb, vb, xp, vp = yb_c, nb, xp_c, np_
    return log

for (mv, Mv, kv, gv) in [(1.0, 3.0, 100.0, 9.8), (0.5, 2.0, 40.0, 9.8)]:
    hv = float(hsol[0].subs({m: mv, M: Mv, k: kv, g: gv}))
    log = simulate(mv, Mv, kv, gv, hv)
    print(f'm={mv} M={Mv} k={kv}: h={hv:.6f}, v0={np.sqrt(2*gv*hv):.6f}')
    for e in log[:4]:
        print('   t=%.6f plate pos=%.2e ball v_in=%.6f plate v_in=%.6f -> ball %.6f plate %.6f' % e)
    per = 2*np.sqrt(2*hv/gv) + np.pi*np.sqrt(Mv/kv)
    print('   period formula', per, ' t(3rd)-t(1st)=', log[2][0]-log[0][0])
