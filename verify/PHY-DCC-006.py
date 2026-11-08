# N31-1: bulb R = R0(1+beta P), R0=30, data 40V->0.50A
import sympy as sp, numpy as np
from scipy.optimize import brentq
R0=30; beta=sp.symbols('beta',positive=True)
b=sp.solve(sp.Eq(80,R0*(1+beta*40*0.5)),beta)[0]; print('beta=',sp.nsimplify(b))
b=sp.Rational(1,12)
I,V=sp.symbols('I V',positive=True)
# R=V/I, P=VI -> V/I = R0(1+beta V I)
Vsol=sp.solve(sp.Eq(V/I,R0*(1+b*V*I)),V)[0]; print('V(I)=',sp.simplify(Vsol))
Im=1/sp.sqrt(b*R0); print('Im=',Im,float(Im))
f=sp.lambdify(I,Vsol)
for i in [0.2,0.4,0.5,0.6]: print(i,f(i))
# bridge: left branch L (top) + R1=50 (bottom), right R2=R3=100; balance R_L=50
RL=50; P=sp.solve(sp.Eq(RL,R0*(1+b*sp.Symbol('P'))))[0]; IL=sp.sqrt(P/RL); print('balance: P=',P,'I=',IL,'E=',IL*(RL+50))
# E such that I_L=0.5 (G removed): E=V_L+50*0.5
E2=f(0.5)+50*0.5; print('E for 0.5A:',E2,' Va(from bottom)=',50*0.5,' Vb=',E2/2,' Vb-Va=',E2/2-25)
# with galvanometer (r_G=10 ohm) connected at E=65: solve nonlinear network numerically, check current direction b->a
def Vbulb(i): return R0*i/(1-b*R0*i*i)
Vbulbf=sp.lambdify(I,Vbulb(I))
def solve(E,rG=10.0):
    # unknowns: Va, Vb (bottom=0, top=E). bulb current from top to a: iL with E-Va=Vbulb(iL)
    def iL_of(Va):
        return brentq(lambda i: float(Vbulbf(i))-(E-Va),0,float(Im)-1e-12) if E-Va>0 else 0.0
    def F(Va):
        # given Va, solve Vb from right node KCL: (E-Vb)/100 = Vb/100 + (Vb-Va)/rG
        Vb=(E/100+Va/rG)/(2/100+1/rG)
        return iL_of(Va)+(Vb-Va)/rG - Va/50, Vb
    Va=brentq(lambda x: F(x)[0],1e-9,E-1e-9); Vb=F(Va)[1]
    return Va,Vb,(Vb-Va)/rG
for E in [40,65,200]: print('E=',E,'Va,Vb,iG(b->a)=',solve(E))
# second bulb replacing R1: midpoint = E/2 always
