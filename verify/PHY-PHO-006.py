# P82-1: 2波長で照らした光電管を電源とみなす．折れ線特性と負荷抵抗
import numpy as np, sympy as sp
from scipy.optimize import brentq, minimize_scalar
hc=1240.0; W=2.0
U1=hc/248-W; U2=hc/496-W; print('U1,U2',U1,U2)
I1,I2=0.40,1.20  # μA
# (1) 一様分布モデルの導出: 法線成分 Kn ~ U(0,Kmax), 到達: Kn>=e u → i=I(1-u/U)
u,U,I=sp.symbols('u U I',positive=True)
print(sp.integrate(I/U,(sp.Symbol('k'),u,U)))
def i(u): return I1*max(0,1-u/U1)+I2*max(0,1-u/U2)
Rstar=U2/i(U2); print('R*',Rstar)
for R in [1.0,3.0]:
    us=brentq(lambda x: x/R-i(x),0,U1); print('R',R,'u',us,'i',us/R)
# 解析式
print('R=1:',1.0*(I1+I2)/(1+1.0*(I1/U1+I2/U2)),' R=3:',3*I1/(1+3*I1/U1))
# 最大電力 (数値で全範囲)
uu=np.linspace(0,U1,300001); P=np.array([x*i(x) for x in uu])
k=np.argmax(P); print('Pmax',P[k],'u',uu[k],'R',uu[k]/i(uu[k]))
a=I1+I2; b=I1/U1+I2/U2; print('seg1 local: u',a/(2*b),'P',a*a/(4*b),'R',(a/(2*b))/(a/2))
print('seg2 local: u',U1/2,'P',I1*U1/4,'R',(U1/2)/(I1/2))
# P(R) の極大を R で数値走査
Rs=np.logspace(-2,2,4001); Pr=[]
for R in Rs:
    us=brentq(lambda x: x/R-i(x),0,U1); Pr.append(us*us/R)
Pr=np.array(Pr); loc=[j for j in range(1,len(Rs)-1) if Pr[j]>Pr[j-1] and Pr[j]>Pr[j+1]]
print('local maxima in R:',[(Rs[j],Pr[j]) for j in loc])

# ---- new (4): switching of the optimal operating point vs I2 ----
I1s,I2s,U1s,U2s=sp.symbols('I1 I2 U1 U2',positive=True)
P1=(I1s+I2s)**2/(4*(I1s/U1s+I2s/U2s)); P2=I1s*U1s/4
print('P2-P1 numerator factor:',sp.factor(sp.simplify((P2-P1)*4*(I1s/U1s+I2s/U2s))))
thr=sp.solve(sp.Eq(P1,P2),I2s); print('thresholds I2:',thr)
r=sp.symbols('r',positive=True)
print('(r-2)-(1-2/r) =',sp.factor((r-2)-(1-2/r)))
R1_thr=sp.simplify((1/(I1s/U1s+I2s/U2s)).subs(I2s,I1s*(U1s/U2s-2)))
print('R1 at threshold * I1/U1 =',sp.factor(sp.simplify(R1_thr*I1s/U1s)))
# numeric brute force: global optimum segment as I2 varies (U1=3, U2=0.5, I1=0.4) -> switch at 1.6
def best(I2v):
    def ii(x): return I1*max(0,1-x/U1)+I2v*max(0,1-x/U2)
    uu=np.linspace(0,U1,300001); P=uu*np.array([ii(x) for x in uu]); k=np.argmax(P)
    return uu[k],P[k],uu[k]/ii(uu[k])
for I2v in [1.2,1.55,1.65,3.0]:
    print(' I2=%.2f -> u*=%.3f P=%.4f R=%.3f'%((I2v,)+best(I2v)))
