# ATM-NRX-101: decay counting
import json, itertools
P = json.load(open('/home/claude/work-G/content/problems/ATM/ATM-NRX-101/problem.json'))
A = {x['id']: x.get('grading') for x in P['answer']['parts']}
sol = [(a, b) for a, b in itertools.product(range(20), range(20)) if 238 - 4*a == 206 and 92 - 2*a + b == 82]
assert sol == [(8, 6)] and A['3'] == '8 [回]' and A['4'] == '6 [回]'
assert (14 + 4 - 17, 7 + 2 - 8) == (1, 1)
assert A['1'] == 'ΔA：-4 ; ΔZ：-2' and A['2'] == 'ΔA：0 ; ΔZ：1'
print('OK')
