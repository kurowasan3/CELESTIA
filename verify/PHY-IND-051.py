# PHY-IND-051 : 3本スポークの車輪の半分だけに磁場（相補性で仕事率一定）
import sympy as sp
import numpy as np
from scipy.integrate import solve_ivp

E, r, s = sp.symbols('E r s', positive=True)

def solve_net(inS):
    """nodal analysis. node O potential 0. spoke k: EMF E (O->Pk, Pk positive) if k in field."""
    V = sp.symbols('V1:4')
    eqs = []
    for k in range(3):
        emf = E if k in inS else 0
        Ispoke = (emf - V[k]) / r          # current O->Pk through spoke
        arcs = sum((V[k] - V[j]) / s for j in range(3) if j != k)  # leaving Pk via arcs
        eqs.append(sp.Eq(Ispoke, arcs))
    sol = sp.solve(eqs, V, dict=True)[0]
    I = [sp.simplify(((E if k in inS else 0) - sol[V[k]]) / r) for k in range(3)]
    P = sp.simplify(sum((E if k in inS else 0) * I[k] for k in range(3)))
    return sol, I, P

for S in [{0}, {0, 1}, {0, 1, 2}]:
    sol, I, P = solve_net(S)
    print(S, {k: sp.simplify(v) for k, v in sol.items()}, I, P)

sol1, I1, P1 = solve_net({0})
assert sp.simplify(I1[0] - 2*E/(3*r+s)) == 0
assert sp.simplify(sol1[sp.Symbol('V1')] - (r+s)*E/(3*r+s)) == 0
assert sp.simplify(sol1[sp.Symbol('V2')] - r*E/(3*r+s)) == 0
sol2, I2, P2 = solve_net({0, 1})
assert sp.simplify(I2[2] + 2*E/(3*r+s)) == 0   # OP3: P3->O で 2E/(3r+s)
assert sp.simplify(I2[0] - E/(3*r+s)) == 0
assert sp.simplify(P1 - P2) == 0 and sp.simplify(P1 - 2*E**2/(3*r+s)) == 0
sol3, I3, P3 = solve_net({0, 1, 2})
assert all(sp.simplify(i) == 0 for i in I3)
# Joule heat check (independent): sum I^2 R
def joule(inS):
    sol, I, P = solve_net(inS)
    V = [sol[sp.Symbol('V%d' % (k+1))] for k in range(3)]
    H = sum(i**2*r for i in I) + sum((V[a]-V[b])**2/s for a, b in [(0,1),(1,2),(2,0)])
    return sp.simplify(H)
assert sp.simplify(joule({0}) - P1) == 0 and sp.simplify(joule({0,1}) - P2) == 0
print('network OK')

# 減速：J dω/dt = -P/ω,  E = B ω a^2 /2 ; 角度依存のスポーク集合を使って数値積分
B, a, rr, ss, J, w0 = 0.8, 0.5, 0.3, 0.7, 2e-3, 30.0
def nin(theta):
    cnt = 0
    for k in range(3):
        ang = (theta + 2*np.pi*k/3) % (2*np.pi)
        if ang > np.pi:   # y<0
            cnt += 1
    return cnt
def rhs(t, y):
    th, w = y
    Ev = 0.5*B*w*a**2
    n = nin(th)
    G11 = 2/(3*rr+ss)
    P = Ev**2*(n*G11 + n*(n-1)*(-G11/2))
    return [w, -P/(J*w) if w > 1e-9 else 0.0]
ev = lambda t, y: y[1]-1e-6
ev.terminal = True
sol = solve_ivp(rhs, [0, 200], [0.1, w0], max_step=1e-3, events=ev, rtol=1e-9, atol=1e-12)
th_num = sol.y[0, -1]-0.1
th_th = 2*J*(3*rr+ss)*w0/(B**2*a**4)
print('theta total num', th_num, 'theory', th_th)
assert abs(th_num-th_th)/th_th < 1e-3
# ω が角度について直線的に減ることを確認
idx = len(sol.t)//2
print('mid: w', sol.y[1, idx], 'pred', w0*(1-(sol.y[0, idx]-0.1)/th_th))
assert abs(sol.y[1, idx]-w0*(1-(sol.y[0, idx]-0.1)/th_th)) < 1e-3
print('ALL OK')
