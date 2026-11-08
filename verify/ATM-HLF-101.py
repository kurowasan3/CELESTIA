# ATM-HLF-101: half-life
import json, math, random
P = json.load(open('/home/claude/work-G/content/problems/ATM/ATM-HLF-101/problem.json'))
A = {x['id']: x.get('grading') for x in P['answer']['parts']}
# (1) memoryless: P(decay in [T,2T] | survived T) = 1/2 ; Monte Carlo check
lam = math.log(2); random.seed(1)
ts = [random.expovariate(lam) for _ in range(400000)]
surv = [t for t in ts if t > 1]
p = sum(1 for t in surv if t <= 2)/len(surv)
assert abs(p - 0.5) < 0.005 and A['1'] == '1/2'
T = 5.7e3
t3 = 3*T
assert float(f'{t3:.1e}') == 1.7e4 and A['3'] == '1.7*10^4'
t5 = T*math.log(5)/math.log(2)
t5_given = T*(1 - 0.30)/0.30
print(t5, t5_given)
assert float(f'{t5:.1e}') == 1.3e4 and float(f'{t5_given:.1e}') == 1.3e4 and A['4'] == '1.3*10^4'
pr = P['predict']
assert 400*0.5**(6/3) == 100 and [c for c in pr['choices'] if c['id'] == pr['answer']][0]['tex'].startswith('$100$')
print('OK')
