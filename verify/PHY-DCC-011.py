# P68-1 N-shaped piecewise-linear element + series R + battery E; load-line, hysteresis
from fractions import Fraction as F
def I_el(V):  # mA, V in volts
    if V<=1: return 10*V
    if V<=3: return 10-4*(V-1)
    return 2+4*(V-3)
def stable_points(E,R):  # R in kOhm -> I(mA)=(E-V)/R
    pts=[]
    # branch1: 10V=(E-V)/R
    V=E/(1+10*R)
    if 0<=V<=1: pts.append(('b1',V,10*V))
    # branch3: 2+4(V-3)=(E-V)/R -> V(4+1/R)=E/R+10
    V=(E/R+10)/(4+1/R)
    if V>=3: pts.append(('b3',V,2+4*(V-3)))
    return pts
R=0.5
print('E=5:',stable_points(5,R))
# sweep up then down with continuity rule
import numpy as np
state='b1'; log=[]
for E in np.concatenate([np.linspace(0,8,80001),np.linspace(8,0,80001)]):
    pts={p[0]:p for p in stable_points(E,R)}
    if state not in pts:
        new=[k for k in pts][0]; log.append((round(E,4),state,'->',new,round(pts[new][2],4))); state=new
print(log)
print('expected up: E=6, I->',14/3,' down: E=4, I->',20/3)
# (4): E=5, from b3 state: raise R until b3 disappears, then back to 0.5
state='b3'; 
for Rv in list(np.linspace(0.5,1.5,10001))+list(np.linspace(1.5,0.5,10001)):
    pts={p[0]:p for p in stable_points(5,Rv)}
    if state not in pts:
        new=list(pts)[0]; print('R=%.4f kOhm: %s->%s'%(Rv,state,new)); state=new
print('final state',state,stable_points(5,0.5))
# check b1 exists for R in [0.5,1.5] at E=5 (needs R>0.4)
print('b1 exists at R=0.5?', any(p[0]=='b1' for p in stable_points(5,0.5)))
# brute-force: intersections via fine grid incl. middle branch, E=5
Vs=np.linspace(0,6,600001); f=np.array([I_el(v) for v in Vs])-(5-Vs)/0.5
print('all crossings E=5:',Vs[np.where(np.diff(np.sign(f)))[0]])
# (5) critical R: E(V)=V+R*I(V) monotone iff 1-4R>0 (R in kOhm)
for Rv in (0.2499,0.25,0.2501):
    Vg=np.linspace(0,6,60001); Eg=Vg+Rv*np.array([I_el(v) for v in Vg])
    print('R=%.4f kOhm: min dE = %.3e'%(Rv,np.min(np.diff(Eg))))
# at R=0.25 sweep up: jump?
state='b1'
for E in np.linspace(0,8,80001):
    pts={p[0]:p for p in stable_points(E,0.25)}
    if state not in pts and pts:
        print('R=0.25 jump at E=%.4f to I=%.4f'%(E,pts[list(pts)[0]][2])); break
