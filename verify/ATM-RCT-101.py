# ATM-RCT-101: alpha decay energy share, activity, power
import json, math, sympy as sp
P = json.load(open('/home/claude/work-G/content/problems/ATM/ATM-RCT-101/problem.json'))
A = {x['id']: x.get('grading') for x in P['answer']['parts']}
m, M, Q, va, vR = sp.symbols('m M Q v_a v_R', positive=True)
s = sp.solve([sp.Eq(m*va, M*vR), sp.Eq(m*va**2/2 + M*vR**2/2, Q)], [va, vR], dict=True)[0]
Ka = sp.simplify(m*s[va]**2/2)
assert sp.simplify(Ka - M*Q/(M+m)) == 0
Qn = 4.78*226/222; print('Q', Qn)
assert round(Qn, 2) == 4.87 and A['2'].endswith('4.87')
K = 4.78*1.6e-13
v = math.sqrt(2*K/6.6e-27); print('v', v)
assert float(f'{v:.1e}') == 1.5e7 and A['3'] == '1.5*10^7 [m/s]'
N = 1.0/226*6.0e23; T = 1.6e3*3.2e7
act = 0.69/T*N; print('activity', act)
assert float(f'{act:.1e}') == 3.6e10 and A['4'] == '3.6*10^10'
Pw = act*Qn*1.6e-13; print('P', Pw)
assert float(f'{Pw:.1e}') == 2.8e-2 and A['5'] == '2.8*10^(-2) [W]'
# rounded-chain consistency
assert float(f'{3.6e10*4.87*1.6e-13:.1e}') == 2.8e-2
print('OK')
