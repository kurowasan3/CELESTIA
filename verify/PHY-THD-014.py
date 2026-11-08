# P32-1 有限熱容量の2物体間で働くカルノー機関
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
n,R,Cv,TH,TL,V1,V2=sp.symbols('n R C_V T_H T_L V1 V2',positive=True)
g1=R/Cv  # gamma-1
V3=V2*(TH/TL)**(1/g1); V4=V1*(TH/TL)**(1/g1)
QH=n*R*TH*sp.log(V2/V1); QL=n*R*TL*sp.log(V3/V4)
print('QL/QH=',sp.simplify(sp.expand_log(QL/QH,force=True)))
# finite bodies: C each, TA0=4T0, TB0=T0. dQA=-C dTA, dQB=C dTB, dQB/dQA=TB/TA
C,T0=1.0,1.0
f=lambda t,y:[-1.0, (y[1]/y[0])]   # parameter: dTA=-ds, dTB = (TB/TA) ds
ev=lambda t,y:y[0]-y[1]; ev.terminal=True
s=solve_ivp(f,[0,10],[4.0,1.0],events=ev,rtol=1e-10,atol=1e-12,dense_output=True)
TA,TB=s.y[0,-1],s.y[1,-1]; print('final TA,TB=',TA,TB)
W=C*(4-TA)-C*(TB-1); print('W/(C T0)=',W,'  eff overall=',W/(C*(4-TA)))
# stop when TA=3
s2=solve_ivp(f,[0,1],[4.0,1.0],rtol=1e-10,atol=1e-12)
print('TA=3 -> TB=',s2.y[1,-1],' W=',(4-3)-(s2.y[1,-1]-1))
# direct contact
print('direct contact final', (4+1)/2)
