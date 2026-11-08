# ATM-CMP-101: Compton. Derive shift (with given approximation) and check exact relativistic too.
import json, math, sympy as sp
P = json.load(open('/home/claude/work-G/content/problems/ATM/ATM-CMP-101/problem.json'))
A = {x['id']: x.get('grading') for x in P['answer']['parts']}
h, m, c, l, lp, phi, v, al = sp.symbols('h m c lambda lambda_p phi v alpha', positive=True)
# momentum: (mv)^2 from components
mv2 = sp.expand((h/l - h/lp*sp.cos(phi))**2 + (h/lp*sp.sin(phi))**2)
mv2 = sp.simplify(sp.trigsimp(mv2))
# energy: m^2 v^2 = 2 m h c (1/l - 1/lp)
eq = sp.Eq(2*m*h*c*(1/l - 1/lp)*l*lp, sp.expand(mv2*l*lp))
lhs = sp.simplify(eq.lhs); rhs = sp.expand(eq.rhs)
# rhs = h^2 (lp/l + l/lp - 2cos) ; replace lp/l + l/lp -> 2
x = sp.symbols('x')
rhs_check = sp.simplify(rhs - h**2*(lp/l + l/lp - 2*sp.cos(phi)))
assert rhs_check == 0
shift = sp.solve(sp.Eq(lhs, h**2*(2 - 2*sp.cos(phi))), lp)[0] - l
assert sp.simplify(shift - h/(m*c)*(1 - sp.cos(phi))) == 0
# numeric
C = 2.43e-12; L = 6.00e-11
lp90 = L + C*(1 - math.cos(math.pi/2))
print('lp90', lp90)
assert float(f'{lp90:.2e}') == 6.24e-11 and A['3'] == '6.24*10^(-11) [m]'
lp180 = L + 2*C
K = 1.24e-6*(1/L - 1/lp180)
print('K eV', K)
assert float(f'{K:.2e}') == 1.55e3 and A['4'] == '1.55*10^3 [eV]'
# exact relativistic Compton gives the same shift; check electron KE nonrel OK (K << 511 keV)
assert K < 0.01*511e3
ratio = 2*C/5.0e-7
print('ratio', ratio)
assert float(f'{ratio:.1e}') == 9.7e-6 and A['5'] == '9.7*10^(-6)'
print('OK')
