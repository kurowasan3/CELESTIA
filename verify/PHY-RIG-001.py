# M04-1: ladder, rough wall (mu1=1/2) and rough floor (mu2); feasibility via LP
import sympy as sp, numpy as np
from scipy.optimize import linprog
L,M,g=1.0,1.0,1.0
c,h=3/5,4/5   # B=(0,0), A=(c,h) at wall x=c
def feasible(mu1,mu2,loads):
    # loads: list of (mass, horizontal coord). unknowns: N,F2 (floor; F2 toward wall +x), R (wall normal, -x), F1 (wall friction, +y can be +/-)
    W=sum(m for m,_ in loads)*g; Mo=sum(m*g*x for m,x in loads)
    # eqs: x: F2 - R =0 ; y: N + F1 - W =0 ; torque about B (ccw +): -R*h ... careful
    # wall force at A=(c,h): (-R, F1): torque = c*F1 - h*(-R) = c F1 + h R ; gravity torque = -Mo
    Aeq=[[0,1,-1,0],[1,0,0,1],[0,0,h,c]]; beq=[0,W,Mo]
    # |F2|<=mu2 N, |F1|<=mu1 R, N,R>=0
    Aub=[[-mu2,1,0,0],[-mu2,-1,0,0],[0,0,-mu1,1],[0,0,-mu1,-1]]; bub=[0]*4
    r=linprog(np.zeros(4),A_ub=Aub,b_ub=bub,A_eq=Aeq,b_eq=beq,bounds=[(0,None),(None,None),(0,None),(None,None)])
    return r.status==0
# (1) min mu2 with ladder only
lo,hi=0,2
for _ in range(60):
    mid=(lo+hi)/2
    if feasible(0.5,mid,[(M,c/2)]): hi=mid
    else: lo=mid
print('(1) mu2_min LP',hi,' pred 6/19=',6/19)
# (2) mu2=1/2, person 2M at distance x from B (horizontal x*c/L)
lo,hi=0,1
for _ in range(60):
    mid=(lo+hi)/2
    if feasible(0.5,0.5,[(M,c/2),(2*M,mid*c)]): lo=mid
    else: hi=mid
print('(2) xmax/L LP',lo,' pred 17/20=',17/20)
# (3) weight m' at B so person reaches A
lo,hi=0,3
for _ in range(60):
    mid=(lo+hi)/2
    if feasible(0.5,0.5,[(M,c/2),(2*M,c),(mid,0)]): hi=mid
    else: lo=mid
print("(3) m' min LP",hi,' pred 9/22=',9/22)
# symbolic concurrency check
m1,m2,x=sp.symbols('mu1 mu2 x',positive=True)
Ls,Ms,gs=sp.symbols('L M g',positive=True)
xp=m2*(sp.Rational(4,5)*Ls+m1*sp.Rational(3,5)*Ls)/(1+m1*m2)
print('x_P',sp.simplify(xp.subs(m1,sp.Rational(1,2))))
print('mu2 min',sp.solve(sp.Eq(xp.subs(m1,sp.Rational(1,2)),sp.Rational(3,10)*Ls),m2))
xG=sp.Rational(3,5)*(Ms*Ls/2+2*Ms*x)/(3*Ms)
print('x max',sp.solve(sp.Eq(xG,xp.subs({m1:sp.Rational(1,2),m2:sp.Rational(1,2)})),x))
mp=sp.symbols('mp',positive=True)
xG3=sp.Rational(3,5)*(Ms*Ls/2+2*Ms*Ls)/(3*Ms+mp)
print("m' min",sp.solve(sp.Eq(xG3,xp.subs({m1:sp.Rational(1,2),m2:sp.Rational(1,2)})),mp))
# (1) indeterminacy: with mu2 = 1/2 and ladder only, show two different admissible solutions
for F1 in [0.0,0.1]:
    # R from torque: c F1 + h R = M g c/2
    R=(M*g*c/2-c*F1)/h; N=M*g-F1; print('F1',F1,'R',R,'N',N,'F2',R,'ok', abs(R)<=N*0.5 and abs(F1)<=0.5*R)
