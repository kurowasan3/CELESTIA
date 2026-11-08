# N25-1 (revised): bridge with capacitors in two arms; A-L: R, L-B: C, A-N: 2C, N-B: 2R, L-N: S + R
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
E,R,C=sp.symbols('E R C',positive=True)
# (1) S open steady: no current; VL=E, VN=0
Q1o,Q2o=C*E,2*C*E
# (2) just after closing S: capacitor voltages continuous
VL0,VN0=E,E-Q2o/(2*C)
print('just after: I_R',(E-VL0)/R,'I_bridge',(VL0-VN0)/R,'I_2R',VN0/(2*R))
iC2N=(VL0-VN0)/R-VN0/(2*R)   # into N plate of C2
Ibat=(E-VL0)/R - iC2N        # battery current out of + terminal = I_R + current into A-plate (= -iC2N)
print('battery current just after',Ibat)
# (3) steady with S closed
I=E/(R+R+2*R); VL=E-R*I; VN=2*R*I; Q1c,Q2c=C*VL,2*C*(E-VN)
print('steady closed: I',I,'Q1',Q1c,'Q2',Q2c)
# (4) open S again
qbat=(Q1o+Q2o)-(Q1c+Q2c); W=E*qbat
U=lambda q1,q2: q1**2/(2*C)+q2**2/(2*2*C)
dU=U(Q1o,Q2o)-U(Q1c,Q2c); print('qbat',qbat,'W',W,'dU',dU,'H',sp.simplify(W-dU))
print('H_R',sp.Rational(1,2)*C*(E-VL)**2,'H_2R',sp.Rational(1,2)*2*C*(E-(E-VN))**2)
# numeric simulation of whole sequence (E=R=C=1)
def make(Sclosed):
    def f(t,y):
        q1,q2,hR,h2R,hb,w=y   # q1: L-plate of C1 ; q2: A-plate of C2
        VL=q1/1.0; VN=1.0-q2/2.0
        iR=(1.0-VL)/1.0; i2R=VN/2.0; ib=(VL-VN)/1.0 if Sclosed else 0.0
        dq1=iR-ib; dq2n=ib-i2R  # N-plate gain
        Ibat=iR-dq2n
        return [dq1,-dq2n,iR**2*1.0,i2R**2*2.0,ib**2*1.0,1.0*Ibat]
    return f
a=solve_ivp(make(True),[0,60],[1,2,0,0,0,0],rtol=1e-11,atol=1e-13,dense_output=True)
y0=a.sol(0); print('num just-after battery current', make(True)(0,y0)[5])
print('num closed final q1,q2',a.y[0,-1],a.y[1,-1])
# find time battery current changes sign
ts=np.linspace(0,10,20001); Ib=[make(True)(0,a.sol(t))[5] for t in ts]
k=next(i for i,v in enumerate(Ib) if v>0); print('battery current sign change near t=',ts[k])
b=solve_ivp(make(False),[0,80],[a.y[0,-1],a.y[1,-1],0,0,0,0],rtol=1e-11,atol=1e-13)
print('num after opening: q1,q2',b.y[0,-1],b.y[1,-1],'H_R',b.y[2,-1],'H_2R',b.y[3,-1],'W',b.y[5,-1],'expected',1/32,1/4,5/4)
