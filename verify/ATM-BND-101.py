# ATM-BND-101: binding energy, D-D fusion
import json
from fractions import Fraction as F
P = json.load(open('/home/claude/work-G/content/problems/ATM/ATM-BND-101/problem.json'))
A = {x['id']: x.get('grading') for x in P['answer']['parts']}
mp, mn, mD, mHe = F('1.0073'), F('1.0087'), F('2.0136'), F('3.0150')
dm = mp + mn - mD
assert dm == F('0.0024') and A['1'] == '0.0024'
B = float(dm*931); print('B', B)
assert round(B, 1) == 2.2 and A['2'].endswith('2.2')
assert round(B/2, 1) == 1.1 and A['3'].endswith('1.1')
Q = float((2*mD - mHe - mn)*931); print('Q', Q)
assert round(Q, 1) == 3.3 and A['4'].endswith('3.3')
# momentum conservation: equal p, K ~ 1/M ; M_He:M_n = 3:1
Kn = Q*3/(3+1); print('Kn', Kn)
assert round(Kn, 1) == 2.4 and A['5'].endswith('2.4')
# with actual masses
Kn2 = Q*float(mHe/(mHe+mn)); assert round(Kn2, 1) == 2.4
print('OK')
