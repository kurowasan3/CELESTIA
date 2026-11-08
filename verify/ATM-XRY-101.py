# ATM-XRY-101: X-ray min wavelength, Bragg
import json, math, sympy as sp
P = json.load(open('/home/claude/work-G/content/problems/ATM/ATM-XRY-101/problem.json'))
A = {x['id']: x.get('grading') for x in P['answer']['parts']}
h, c, e, V, lam = sp.symbols('h c e V lambda', positive=True)
lmin = sp.solve(sp.Eq(e*V, h*c/lam), lam)[0]
assert sp.simplify(lmin - h*c/(e*V)) == 0
v = float(lmin.subs({h: 6.6e-34, c: 3.0e8, e: 1.6e-19, V: 3.0e4}))
print('lmin', v)
assert float(f'{v:.1e}') == 4.1e-11 and A['2'] == '4.1*10^(-11) [m]'
s1 = 1.5e-10 / (2*2.5e-10)
assert abs(s1 - 0.30) < 1e-12 and A['3'] == '0.30'
nmax = max(n for n in range(1, 20) if n*1.5e-10/(2*2.5e-10) <= 1)
assert nmax == 3 and A['4'] == '3'
r = float(lmin.subs(V, 6.0e4) / lmin.subs(V, 3.0e4))
assert abs(r - 0.5) < 1e-12 and A['5'] == '(1/2)倍'
assert 1.5e-10 > v  # characteristic line longer than cutoff
print('OK')
