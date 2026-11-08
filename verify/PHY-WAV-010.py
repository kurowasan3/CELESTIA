# T20073-1 elementary-wave (phasor) sum for openings in a wall, far point P in direction theta
import numpy as np
import sympy as sp

w0 = 1.0                     # w0 = lambda/sin(theta)
kap = 2*np.pi/w0             # phase per unit length along the wall


def amp(intervals, N=200001):
    tot = 0j
    for lo, hi in intervals:
        if hi <= lo:
            continue
        x = np.linspace(lo, hi, N)
        tot += np.trapezoid(np.exp(1j*kap*x), x)
    return abs(tot)


A1 = amp([(0, w0/2)])
# (1),(2): A(w)/A1 = |sin(pi w/w0)|
for w in [0.1, 0.25, 0.5, 0.8, 1.0, 1.3]:
    print('w=%.2f  sim %.6f  claim %.6f' % (w, amp([(0, w)])/A1, abs(np.sin(np.pi*w/w0))))
# smallest w>0 with zero amplitude
ws = np.linspace(0.01, 1.5, 1491); As = [amp([(0, w)], 4001) for w in ws]
print('first zero near', ws[int(np.argmin(As[:1200]))])
# (3): opening [0,w0], block [p,p+h]
for p, h in [(0.0, 0.3), (0.2, 0.3), (0.55, 0.3), (0.1, 0.5), (0.3, 0.7)]:
    print('p=%.2f h=%.2f sim %.6f claim %.6f' % (p, h, amp([(0, p), (p + h, w0)])/A1, abs(np.sin(np.pi*h/w0))))
# (4): opening [0, 3w0/2], block [p, p+w0/2], 0<=p<=w0
for p in [0, 0.1, 0.25, 0.5, 0.7, 1.0]:
    print('p=%.2f sim %.6f claim 2|sin(pi p/w0)| = %.6f' % (p, amp([(0, p), (p + 0.5, 1.5)])/A1, 2*abs(np.sin(np.pi*p/w0))))
print('unblocked [0,3w0/2]:', amp([(0, 1.5)])/A1)
# symbolic check of the chord formula
x, w, k = sp.symbols('x w k', positive=True)
S = sp.integrate(sp.exp(sp.I*k*x), (x, 0, w))
print('|S|^2 simplify:', sp.simplify(sp.expand_complex(S*sp.conjugate(S)) - (2*sp.sin(k*w/2)/k)**2))
# (5) brute force: any single strip [p,p+h] inside [0,1.5] -> max amplitude
best = (0, None)
for p in np.linspace(0, 1.5, 61):
    for h in np.linspace(0, 1.5 - p, 61):
        A = amp([(0, p), (p + h, 1.5)], 801)/A1
        if A > best[0]:
            best = (A, (round(p, 3), round(h, 3)))
print('brute-force max over single blocked strip:', best)
