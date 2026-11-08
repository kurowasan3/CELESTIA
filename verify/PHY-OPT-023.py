# K20233-1: chain of lossy polarizers rotating polarization by Theta
import numpy as np, sympy as sp
N,Th,t=sp.symbols('N Theta t',positive=True)
# Jones-type numeric simulation
def chain(Nn,tt,Thn):
    E=np.array([0.0,1.0]) # along y (components x,y)
    I0=1.0
    for k in range(1,Nn+1):
        ang=k*Thn/Nn
        u=np.array([np.sin(ang),np.cos(ang)])   # axis at angle from y
        E=np.sqrt(tt)*np.dot(E,u)*u
    return np.dot(E,E)
Thn=np.pi/2; tt=0.98
vals=[(Nn,chain(Nn,tt,Thn)) for Nn in range(1,40)]
best=max(vals,key=lambda p:p[1]); print('best N',best, ' neighbors',vals[best[0]-2],vals[best[0]])
for Nn in [3,10,11,12]:
    print(Nn,chain(Nn,tt,Thn), tt**Nn*np.cos(Thn/Nn)**(2*Nn))
Nstar=Thn/np.sqrt(np.log(1/tt)); print('N* approx',Nstar,' Imax approx',np.exp(-2*Thn*np.sqrt(np.log(1/tt))))
# series check of cos^(2N)(Theta/N)
x=sp.symbols('x',positive=True)  # x=1/N
expr=sp.log(sp.cos(Th*x))*2/x
print(sp.series(expr,x,0,3))
# t=1 lossless, N=2 polarizers each 45deg
print('ideal N=2',chain(2,1.0,np.pi/2),' N=3',chain(3,1.0,np.pi/2), 27/64)
# also: Theta=pi/2, t=0.98 single crossing impossible: N=1 -> 0
