# M46-1: heater calorimetry with two runs, then mixing with phase change
import sympy as sp
P,m,th0,k1,k2,tau,cw,K,th,t=sp.symbols('P m theta0 k1 k2 tau c_w K theta t',positive=True)
ci,C,L=sp.symbols('c_i C L',positive=True)
sol=sp.solve([sp.Eq(P,(m*ci+C)*k1),sp.Eq(P,(2*m*ci+C)*k2),sp.Eq(P*tau,m*L)],[ci,C,L],dict=True)[0]
print('ci',sp.simplify(sol[ci]),' C',sp.simplify(sol[C]),' L',sol[L])
assert sp.simplify(sol[ci]-P*(1/k2-1/k1)/m)==0
assert sp.simplify(sol[C]-P*(2/k1-1/k2))==0
H0=(m*sol[ci]+sol[C])*th0   # heat to bring ice+container to 0
print('H0',sp.simplify(H0))
assert sp.simplify(H0-P*th0/k1)==0
Hmelt=m*sol[L]
lo=sp.simplify(H0/K); hi=sp.simplify((H0+Hmelt)/K)
print('range',lo,hi)
# melted mass
x=sp.simplify((K*th-H0)/sol[L]); print('melted',x)
# final temp above 0
tf=sp.solve(sp.Eq(K*(th-t),H0+Hmelt+(m*cw+sol[C])*t),t)[0]
print('tf',sp.simplify(tf))
assert sp.simplify(tf-(K*th-P*th0/k1-P*tau)/(K+m*cw+P*(2/k1-1/k2)))==0
# numeric sanity: P=500W, m=100g, ci=2.1, C=50 J/K, L=334 J/g
vals={P:500,m:100}
ci_n,C_n,L_n=2.1,50,334
k1n=500/(100*ci_n+C_n); k2n=500/(200*ci_n+C_n); taun=100*L_n/500
print('k1,k2,tau',k1n,k2n,taun, 'k1/2<k2<k1', k1n/2<k2n<k1n)
r=sol[ci].subs({P:500,m:100,k1:k1n,k2:k2n}); print('ci back',float(r))
r=sol[C].subs({P:500,m:100,k1:k1n,k2:k2n}); print('C back',float(r))
# brute-force equilibrium for a block: K=400 J/K, theta=90, theta0=10, cw=4.2
def final(Kb,thb,th0n=10,cwn=4.2):
    # energy function: enthalpy of ice/water+container relative to ice at 0
    import numpy as np
    from scipy.optimize import brentq
    def Hsys(T):  # T in C, heat content of sample+container relative to ice at -th0
        if T<0: return (100*ci_n+C_n)*(T+th0n)
        return (100*ci_n+C_n)*th0n + 100*L_n + (100*cwn+C_n)*T
    f=lambda T: Hsys(T)-((100*ci_n+C_n)*th0n*0) - Kb*(thb-T)
    # check plateau
    if Hsys(0)-0 <= Kb*thb <= Hsys(0)+0 and False: pass
    H0n=(100*ci_n+C_n)*th0n
    if Kb*thb<H0n: return brentq(lambda T:(100*ci_n+C_n)*(T+th0n)-Kb*(thb-T),-th0n,0),None
    if Kb*thb<=H0n+100*L_n: return 0,(Kb*thb-H0n)/L_n
    return brentq(lambda T:H0n+100*L_n+(100*cwn+C_n)*T-Kb*(thb-T),0,thb),100
for Kb,thb in [(400,10),(400,50),(400,95),(800,90)]:
    num=final(Kb,thb)
    sub={P:500,m:100,k1:k1n,k2:k2n,tau:taun,th0:10,cw:4.2,K:Kb,th:thb}
    lo_n=float(lo.subs(sub)); hi_n=float(hi.subs(sub))
    if thb>hi_n: print(Kb,thb,'num',num,'formula tf',float(tf.subs(sub)))
    elif thb>=lo_n: print(Kb,thb,'num',num,'formula melted',float(x.subs(sub)))
    else: print(Kb,thb,'num',num,'below range')
print('OK')
