import sympy as sp
C,C0,CY,Q,q0=sp.symbols('C C_0 C_Y Q q_0',positive=True)
VE,VY=sp.symbols('V_E V_Y')
# nodes: E (charge -q0), Y (charge Q); G ground
eqs=[sp.Eq(C*(VE-VY)+C0*VE,-q0), sp.Eq(C*(VY-VE)+CY*VY,Q)]
s=sp.solve(eqs,[VE,VY])
qL=sp.simplify(C0*s[VE]); print('q_L(C)=',sp.factor(qL))
print('C at close:',sp.solve(sp.Eq(qL,0),C))
print('C->oo:',sp.limit(qL,C,sp.oo))
print('C=0:',qL.subs(C,0))
print('dqL/dC:',sp.factor(sp.diff(qL,C)))
# grounding at C1
C1=sp.symbols('C_1',positive=True)
VYg=Q/(C1+CY); qE=C1*(0-VYg)
print('E charge after grounding:',sp.simplify(qE))
# flow through S: change of E charge
print('flow into E:',sp.factor(qE-(-q0)))
print('zero flow at:',sp.solve(sp.Eq(qE,-q0),C1))
# then Y removed: leaf charge = qE (all on L side, C->0)
import random
ok=True
for _ in range(2000):
    Qv=random.uniform(1,5); q0v=random.uniform(0.05,0.95)*Qv; CYv=random.uniform(0.1,3); C1v=random.uniform(0,20)
    Cs=CYv*q0v/(Qv-q0v)
    if abs(C1v-Cs)<1e-6: continue
    ok &= ((C1v*Qv/(C1v+CYv)>q0v)==(C1v>Cs))
print('final |qL|>q0 iff C1>C* :',ok)
# Case A of book check: initial neutral, CY->0 limit gives q_L->Q
sA=sp.solve([sp.Eq(C*(VE-VY)+C0*VE,0), sp.Eq(C*(VY-VE)+CY*VY,Q)],[VE,VY])
print('neutral case qL:',sp.factor(C0*sA[VE]))
