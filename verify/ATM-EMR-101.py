# ATM-EMR-101: e/m from magnetic deflection
import json, sympy as sp
P = json.load(open('/home/claude/work-G/content/problems/ATM/ATM-EMR-101/problem.json'))
A = {x['id']: x for x in P['answer']['parts']}
m, e, V, B, r, v = sp.symbols('m e V B r v', positive=True)
vs = sp.solve(sp.Eq(m*v**2/2, e*V), v)[0]
assert sp.simplify(vs - sp.sqrt(2*e*V/m)) == 0
rs = sp.solve(sp.Eq(m*v**2/r, e*v*B), r)[0]
assert sp.simplify(rs - m*v/(e*B)) == 0
# eliminate v
em = sp.solve(sp.Eq(rs.subs(v, vs), r), e)[0] / m
assert sp.simplify(em - 2*V/(B**2*r**2)) == 0
val = 2*2.0e2/((1.0e-3)**2*(4.8e-2)**2)
print('e/m =', val)
assert float(f'{val:.1e}') == 1.7e11 and A['4']['grading'] == '1.7*10^11'
print('OK')
