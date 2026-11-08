# T20112-1: Dickson charge pump with diodes (switch + resistor model, optional forward drop VD)
import numpy as np
from scipy.integrate import solve_ivp
def settle(q,E,X,Y,N,VD,C=1.0,r=1e-3,T=0.5):
    # nodes 1..N pump caps (odd->X, even->Y), node N+1 output cap to G
    def volts(q):
        V=np.zeros(N+2); V[0]=E
        for k in range(1,N+1):
            V[k]=q[k-1]/C + (X if k%2==1 else Y)
        V[N+1]=q[N]/C
        return V
    def rhs(t,q):
        V=volts(q); i=[max(V[k]-V[k+1]-VD,0)/r for k in range(N+1)]
        return [i[k]-i[k+1] for k in range(N)]+[i[N]]
    s=solve_ivp(rhs,[0,T],q,method='LSODA',rtol=1e-10,atol=1e-12)
    qf=s.y[:,-1]; return qf, volts(qf)
def run(N,E,V0,VD,cycles):
    q=np.zeros(N+1)
    q,V=settle(q,E,0,0,N,VD); hist=[('init',V.copy())]
    for c in range(cycles):
        q,V=settle(q,E,V0,0,N,VD); hist.append(('a',V.copy()))
        q,V=settle(q,E,0,V0,N,VD); hist.append(('b',V.copy()))
    return hist
E,V0=5.0,6.0
h=run(2,E,V0,0.0,1)
for tag,V in h: print(tag,np.round(V,4))
print('expected init all',E,' after a: E+V0/3=',E+V0/3,' after b: out E+5V0/6=',E+5*V0/6,' node1=',E)
for N in (2,3,5):
    for VD in (0.0,0.7):
        hh=run(N,5.0,5.0,VD,400)
        print('N',N,'VD',VD,'out',round(hh[-1][1][N+1],4),'formula',5+N*5-(N+1)*VD)
print('N needed for >=30V:',[n for n in range(1,12) if 5+5*n-0.7*(n+1)>=30][0], 'out N=6:',5+30-0.7*7, 'N=5:',5+25-0.7*6)
