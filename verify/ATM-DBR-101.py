# ATM-DBR-101: electron diffraction
import json, math, sympy as sp
P = json.load(open('/home/claude/work-G/content/problems/ATM/ATM-DBR-101/problem.json'))
A = {x['id']: x.get('grading') for x in P['answer']['parts']}
h, m, e, V, c = sp.symbols('h m e V c', positive=True)
lam = h/sp.sqrt(2*m*(e*V))   # p = sqrt(2mK)
num = dict(h=6.6e-34, m=9.1e-31, e=1.6e-19, c=3.0e8)
l150 = float(lam.subs({h: num['h'], m: num['m'], e: num['e'], V: 150}))
print('lambda', l150)
assert float(f'{l150:.1e}') == 1.0e-10 and A['2'] == '1.0*10^(-10) [m]'
d = 1.0e-10/0.46     # rounded lambda as in solution
d_full = l150/0.46
print('d', d, d_full)
assert float(f'{d:.1e}') == 2.2e-10 and float(f'{d_full:.1e}') == 2.2e-10 and A['3'] == 'd：2.2*10^(-10) [m]'
s600 = float(lam.subs({h: num['h'], m: num['m'], e: num['e'], V: 600}))/d_full
print('sin600', s600)
assert abs(s600 - 0.23) < 0.005 and A['4'] == '0.23'
ratio = sp.simplify((h*c/(e*V))/lam)
assert sp.simplify(ratio - sp.sqrt(2*m*c**2/(e*V))) == 0
r = float(ratio.subs({h: num['h'], m: num['m'], e: num['e'], c: num['c'], V: 150}))
print('ratio', r)
assert round(r) == 83 and A['6'] == '83倍'
print('OK')
