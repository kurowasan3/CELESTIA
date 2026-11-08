# K20222-1 verification: long finite ladder (each rail segment r=2R/3, rungs R), emf E in one rung
# nodal analysis with ideal LED modeled as stiff voltage source Vd (checked conducting: current>0)
import numpy as np, sympy as sp
Rs=sp.symbols('R_inf',positive=True)
print('R_inf (R=1):',sp.solve(sp.Eq(Rs,sp.Rational(4,3)+Rs/(1+Rs)),Rs))
R=1.0; r=2*R/3; Nc=80; n=2*Nc+1; idx=lambda k:k+Nc
def solve(emf, led=None, Vd=0.0, gL=1e9):
    G=np.zeros((2*n,2*n)); I=np.zeros(2*n)
    def add(a,b,res):
        G[a,a]+=1/res; G[b,b]+=1/res; G[a,b]-=1/res; G[b,a]-=1/res
    for k in range(-Nc,Nc):
        add(idx(k),idx(k+1),r); add(n+idx(k),n+idx(k+1),r)
    for k in range(-Nc,Nc+1):
        e=emf.get(k,0.0); t=idx(k); b=n+idx(k); add(t,b,R); I[t]+=e/R; I[b]-=e/R
    add(idx(-Nc),n+idx(-Nc),2*R); add(idx(Nc),n+idx(Nc),2*R)   # terminate with R_inf
    if led is not None:
        t=idx(led); b=n+idx(led); add(t,b,1/gL); I[t]+=gL*Vd; I[b]-=gL*Vd
    g0=n+idx(0); G[g0,:]=0; G[g0,g0]=1; I[g0]=0
    x=np.linalg.solve(G,I); T=x[:n]; B=x[n:]
    V=lambda k:T[idx(k)]-B[idx(k)]
    up=lambda k:(emf.get(k,0.0)-V(k))/R          # rung current bottom->top
    iled=(lambda: gL*(Vd-V(led))) if led is not None else None  # current through LED top->bottom
    return V,up,iled
E=1.0
V,up,_=solve({0:E})
print('I0=',up(0),'(E/2R=0.5)  V_n=',[round(V(k),6) for k in range(5)],' (E/2)3^-n=',[round(0.5/3**k,6) for k in range(5)])
print('V_-n symmetric:',[round(V(-k),6) for k in range(3)])
# LED at rung 0 (same as emf rung) and at rung 1, Vd=0.1E
for led in (0,1,2):
    Vd=0.1
    V2,up2,il=solve({0:E},led=led,Vd=Vd)
    pred=(0.5/3**led - Vd)/(R/2)
    print('LED at rung',led,': LED current',il(),' Thevenin pred',pred,' I0=',up2(0))
# power with LED lit at emf rung: P = E*I0, pred E(E-Vd)/R
V3,up3,il3=solve({0:E},led=0,Vd=0.1); print('P=',E*up3(0),' pred',E*(E-0.1)/R, ' LED I pred (E-2Vd)/R=',(E-0.2)/R, il3())
