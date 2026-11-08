# N47-1: resonant charging with diode, then LC transfer C1 -> C2=kC through L2 and diode
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
E,L1,L2,C,k,t=sp.symbols('E L1 L2 C k t',positive=True)
q=sp.Function('q')
# stage 1: L1 q'' + q/C = E
sol=sp.dsolve(L1*q(t).diff(t,2)+q(t)/C-E,q(t),ics={q(0):0,q(t).diff(t).subs(t,0):0}).rhs
print('stage1 q(t)=',sp.simplify(sol))
i1=sp.diff(sol,t); T1=sp.pi*sp.sqrt(L1*C)
print('i at T1:',sp.simplify(i1.subs(t,T1)),' V at T1:',sp.simplify(sol.subs(t,T1)/C),' imax:',sp.simplify(i1.subs(t,T1/2)))
print('battery work - stored:',sp.simplify(E*sol.subs(t,T1)-(sol.subs(t,T1))**2/(2*C)))
# stage 2
Q0=2*C*E
sol2=sp.dsolve(L2*q(t).diff(t,2)+(q(t)-Q0)/C+q(t)/(k*C),q(t),ics={q(0):0,q(t).diff(t).subs(t,0):0}).rhs
sol2=sp.simplify(sol2); print('stage2 q(t)=',sol2)
T2=sp.pi*sp.sqrt(L2*C*k/(1+k))
print('i at T2:',sp.simplify(sp.diff(sol2,t).subs(t,T2)))
qf=sp.simplify(sol2.subs(t,T2))
V1=sp.simplify((Q0-qf)/C); V2=sp.simplify(qf/(k*C))
print('V1=',sp.factor(V1),' V2=',sp.factor(V2))
print('fraction to C2:',sp.factor(sp.simplify((k*C*V2**2/2)/(2*C*E**2))))
print('energy check:',sp.simplify(C*V1**2/2+k*C*V2**2/2-2*C*E**2))
print('k=3:',V1.subs(k,3),V2.subs(k,3))
print('resistor case: common V',sp.simplify(Q0/(C+k*C)),' loss fraction',sp.simplify((2*C*E**2-(C+k*C)*(Q0/(C+k*C))**2/2)/(2*C*E**2)))
# stage-2 current max: (q_e)*omega2
imax2=sp.simplify(sp.diff(sol2,t).subs(t,T2/2)); print('imax2',imax2)
# numeric with ideal diode: integrate until current returns to 0
p=dict(E=5.0,L1=0.2,L2=0.1,C=1e-3,k=3.0)
def f1(t,y):
    qq,ii=y; return [ii,(p['E']-qq/p['C'])/p['L1']]
ev=lambda t,y:y[1]; ev.terminal=True; ev.direction=-1
s=solve_ivp(f1,[0,1],[0,0],events=ev,rtol=1e-11,atol=1e-14,max_step=1e-4)
Vc1=s.y_events[0][0][0]/p['C']; print('num stage1: t',s.t_events[0][0],'V',Vc1, 'imax',s.y[1].max())
def f2(t,y):
    qq,ii=y; return [ii,((2*p['C']*p['E']-qq)/p['C']-qq/(p['k']*p['C']))/p['L2']]
s=solve_ivp(f2,[0,1],[0,0],events=ev,rtol=1e-11,atol=1e-14,max_step=1e-4)
qf_n=s.y_events[0][0][0]
print('num stage2: t',s.t_events[0][0],'V1',(2*p['C']*p['E']-qf_n)/p['C'],'V2',qf_n/(p['k']*p['C']))
print('formula: t',float(T2.subs({L2:p['L2'],C:p['C'],k:p['k']})),'V1',float(V1.subs({E:p['E'],k:p['k']})),'V2',float(V2.subs({E:p['E'],k:p['k']})))
