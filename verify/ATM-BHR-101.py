# ATM-BHR-101: Bohr model, hydrogen-like ions
import json, sympy as sp
P = json.load(open('/home/claude/work-G/content/problems/ATM/ATM-BHR-101/problem.json'))
A = {x['id']: x.get('grading') for x in P['answer']['parts']}
m, e, k0, h, Z, r, v = sp.symbols('m e k_0 h Z r v', positive=True)
n = sp.symbols('n', positive=True, integer=True)
sol = sp.solve([sp.Eq(m*v**2/r, k0*Z*e**2/r**2), sp.Eq(2*sp.pi*r, n*h/(m*v))], [r, v], dict=True)[0]
rn = sol[r]; vn = sol[v]
assert sp.simplify(rn - n**2*h**2/(4*sp.pi**2*k0*m*Z*e**2)) == 0
En = sp.simplify(m*vn**2/2 - k0*Z*e**2/rn)
assert sp.simplify(En + 2*sp.pi**2*k0**2*m*Z**2*e**4/(n**2*h**2)) == 0
# numeric check of 13.6 eV with SI constants
E1 = float((2*sp.pi**2*k0**2*m*e**4/h**2).subs({k0: 8.988e9, m: 9.109e-31, e: 1.602e-19, h: 6.626e-34}))/1.602e-19
assert abs(E1 - 13.6) < 0.05
dE = 13.6*(sp.Rational(1, 4) - sp.Rational(1, 9))
lam = 1.24e-6/float(dE); print('lambda', lam)
assert float(f'{lam:.2e}') == 6.56e-7 and A['3'] == '6.56*10^(-7) [m]'
# He+: Z=2, find upper level when lower is 4
up = [nu for nu in range(5, 50) if abs(4*13.6*(1/16 - 1/nu**2) - float(dE)) < 1e-9]
assert up == [6] and A['4'] == '6'
assert abs(4*13.6 - 54.4) < 1e-12 and A['5'] == '54.4 [eV]'
print('OK')
