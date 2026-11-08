# N24-1 verification: capacitor network with two batteries, energy bookkeeping
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
C,E=sp.symbols('C E',positive=True)
x=sp.symbols('x')
# Step1: S1 closed (P=3E), Q node isolated (charge 0)
VQ1=sp.solve(C*(x-3*E)+2*C*x,x)[0]; print('VQ1',VQ1)
Qa1=C*(3*E-VQ1); Qb1=2*C*VQ1; U1=sp.Rational(1,2)*C*(3*E-VQ1)**2+sp.Rational(1,2)*2*C*VQ1**2
print('Qa1,Qb1,U1',Qa1,Qb1,U1)
# Procedure A: S1,S2 closed: P=3E, Q=2E
Qa2=C*(3*E-2*E); Qb2=2*C*2*E; U2=sp.Rational(1,2)*C*E**2+sp.Rational(1,2)*2*C*(2*E)**2
qS2=(-Qa2+Qb2)-(-Qa1+Qb1); qS1=Qa2-Qa1
W2=2*E*qS2; W1=3*E*qS1; H=W1+W2-(U2-U1)
print('A: qS2',qS2,'W1',W1,'W2',W2,'dU',U2-U1,'H',sp.simplify(H))
# then open S1,S2, close S3 (P=0), Q isolated charge -Qa2+Qb2
qQ=-Qa2+Qb2
VQ3=sp.solve(C*(x-0)+2*C*x-qQ,x)[0]; U3=sp.Rational(1,2)*C*VQ3**2+sp.Rational(1,2)*2*C*VQ3**2
print('VQ3',VQ3,'H3',sp.simplify(U2-U3),'charge through S3',(C*(0-VQ3))-Qa2)
# Procedure B: from step1 open S1 (P isolated, plate charge Qa1), close S2 (Q=2E)
VP=sp.solve(C*(x-2*E)-Qa1,x)[0]
QbB=2*C*2*E; QaB=C*(VP-2*E); qS2B=(-QaB+QbB)-(-Qa1+Qb1)
UB=sp.Rational(1,2)*C*(VP-2*E)**2+sp.Rational(1,2)*2*C*(2*E)**2
HB=2*E*qS2B-(UB-U1); print('B: VP',VP,'qS2',qS2B,'H',sp.simplify(HB))
# numerical ODE check of procedure A heat (C=1,E=1,R1=1,R2=2.3)
Cn,En,R1,R2,R3=1.0,1.0,1.0,2.3,0.7
def rhsA(t,y):
    VP,VQ,h=y  # node voltages; charges: Ca=C(VP-VQ), Cb=2C VQ
    i1=(3*En-VP)/R1; i2=(2*En-VQ)/R2
    # node P: i1 = C d(VP-VQ)/dt ; node Q: i2 + C d(VP-VQ)/dt = 2C dVQ/dt
    dVQ=(i2+i1)/(2*Cn); dVP=i1/Cn+dVQ
    return [dVP,dVQ,R1*i1**2+R2*i2**2]
s=solve_ivp(rhsA,[0,80],[3*En*2/3*0+3*En- (1.0),1.0,0],rtol=1e-10,atol=1e-12)
# initial: VQ1=E, VP=3E
s=solve_ivp(rhsA,[0,80],[3.0,1.0,0],rtol=1e-10,atol=1e-12)
print('numeric A heat',s.y[2,-1],'final',s.y[0,-1],s.y[1,-1])
def rhs3(t,y):
    VP,VQ,h=y  # S3: P to ground via R3, Q isolated: charge conservation
    i3=-VP/R3  # current into P from ground
    # Q isolated: C(VQ-VP)+2C VQ const -> 3dVQ = dVP
    # node P: i3 = C d(VP-VQ)/dt = C(dVP - dVP/3) = (2/3)C dVP
    dVP=i3/(2/3*Cn); dVQ=dVP/3
    return [dVP,dVQ,R3*i3**2]
s3=solve_ivp(rhs3,[0,80],[3.0,2.0,0],rtol=1e-10,atol=1e-12)
print('numeric S3 heat',s3.y[2,-1],'VQ',s3.y[1,-1])
def rhsB(t,y):
    VP,VQ,h=y  # P isolated: d(VP-VQ)=0
    i2=(2*En-VQ)/R2; dVQ=i2/(2*Cn); return [dVQ,dVQ,R2*i2**2]
sB=solve_ivp(rhsB,[0,80],[3.0,1.0,0],rtol=1e-10,atol=1e-12)
print('numeric B heat',sB.y[2,-1],'VP',sB.y[0,-1])
