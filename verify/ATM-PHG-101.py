# ATM-PHG-101: K0-f graph
import json
P = json.load(open('/home/claude/work-G/content/problems/ATM/ATM-PHG-101/problem.json'))
A = {x['id']: x['grading'] for x in P['answer']['parts']}
f1, K1, f2, K2 = 5.0e14, 0.0, 10.0e14, 2.0
slope_eV = (K2 - K1) / (f2 - f1)          # eV s
h = slope_eV * 1.6e-19
f0 = f1 - K1 / slope_eV
W = slope_eV * f0
fB = 3.0 / slope_eV
print(f0, h, W, fB)
assert abs(f0 - 5.0e14) < 1 and A['1'] == '5.0*10^14 [Hz]'
assert float(f'{h:.1e}') == 6.4e-34 and A['2'] == '6.4*10^(-34)'
assert abs(W - 2.0) < 1e-12 and A['3'] == '2.0 [eV]'
assert abs(fB - 7.5e14) < 1 and A['4'].endswith('7.5*10^14 [Hz]')
# graph figure: points used in fig match the data (scale 20px/1e14Hz, 40px/eV, origin (60,200))
fig = json.load(open('/home/claude/work-G/content/problems/ATM/ATM-PHG-101/q-graph.fig.json'))
dots = [it['at'] for it in fig['items'] if it['t'] == 'dot']
assert [60 + 20*5, 200] in dots and [60 + 20*10, 200 - 40*2] in dots
line = [it for it in fig['items'] if it['t'] == 'line' and it.get('style') == 'thick'][0]
(x0, y0), (x1, y1) = line['from'], line['to']
assert abs((y0 - y1) / (x1 - x0) - (40*2)/(20*5)) < 1e-9
print('OK')
