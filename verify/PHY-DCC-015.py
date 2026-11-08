# Q158-1 有限はしご形回路: 末端からの逆算(フィボナッチ), 誤差の漸化式, 整合終端
import sympy as sp, numpy as np
R,X,I0=sp.symbols('R X I0',positive=True)
def ladder_input(n,Xend):
    # n sections: series R then shunt; last shunt replaced by Xend
    Z=Xend
    for i in range(n):
        Z = R + (R*Z/(R+Z) if i>0 else Z)
    return sp.simplify(Z)
# all-R ladder
for n in range(1,7): print(n, sp.nsimplify(ladder_input(n,R)/R))
# back-calculation 5 sections, last rung current I0
v=I0*R; i=I0
for sec in range(5):
    if sec>0: i=i+v/R   # add rung current at this node
    v=v+i*R            # series resistor
print('5 sections: V=',sp.simplify(v),' I=',sp.simplify(i),' Rin=',sp.simplify(v/i))
# numerical check by solving Kirchhoff with numpy for 5 sections
import itertools
def solve_ladder(n,Rv=1.0,Xv=1.0,Vin=1.0):
    # nodes 1..n ; node0 = input at Vin ; ground 0V
    A=np.zeros((n,n)); b=np.zeros(n)
    for k in range(n):
        g_sh = 1/Xv if k==n-1 else 1/Rv
        A[k,k]+=1/Rv+g_sh; 
        if k==0: b[k]+=Vin/Rv
        else: A[k,k-1]-=1/Rv; A[k-1,k]-=1/Rv; A[k-1,k-1]+=1/Rv
    V=np.linalg.solve(A,b)
    Iin=(Vin-V[0])/Rv; Ilast=V[-1]/Xv
    return Vin/Iin, Vin/Ilast
print('numeric 5 sections: Rin, V/I_last =',solve_ladder(5))
phi=(1+np.sqrt(5))/2
Rn=[float(ladder_input(n,R).subs(R,1)) for n in range(1,9)]
print('rel err:',[ (r-phi)/phi for r in Rn])
# exact error relation
Rn_=sp.symbols('Rn',positive=True); ph=(1+sp.sqrt(5))/2
lhs=(R+R*Rn_/(R+Rn_))-ph*R; rhs=R**2*(Rn_-ph*R)/((R+Rn_)*(R+ph*R))
print('error relation:',sp.simplify(lhs-rhs))
# matched termination
print('X=0:',[sp.nsimplify(ladder_input(n,sp.Integer(0))/R) for n in range(1,6)])
# matched termination: the last shunt replaced by X; independence of n requires X = R/phi
Xs=sp.symbols('Xs',positive=True)
print('solve Z1=Z2:',sp.solve(sp.Eq(ladder_input(1,Xs),ladder_input(2,Xs)),Xs))
for n in range(1,6): print('X=R/phi, n=',n, sp.nsimplify(sp.simplify(ladder_input(n,R*(sp.sqrt(5)-1)/2)/R)))
print('numeric check X=R/phi n=3:',solve_ladder(3,1.0,(np.sqrt(5)-1)/2)[0])
