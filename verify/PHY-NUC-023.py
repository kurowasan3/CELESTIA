# S0610-1 accretion vs fusion energy, species-dependent ejection
import sympy as sp
G,M,m,R,Q,x,N=sp.symbols('G M m R Q x N',positive=True)
K=G*M*m/R
eta=sp.simplify(N*K/(N/4*Q)); print("eta",eta)
# energy per atom after burning fraction x of N H-atoms
nH=N*(1-x); nHe=N*x/4; E=N*x/4*Q
e=sp.simplify(E/(nH+nHe)); print("e(x)",e)
k=sp.symbols('k',positive=True)
xH=sp.solve(sp.Eq(e,k),x); print("x for e=k:",xH)
print("de/dx",sp.simplify(sp.diff(e,x)))
num=dict(G=6.67e-11,m=1.67e-27,Q=4.28e-12)
def vals(Mv,Rv):
    Kv=num['G']*Mv*num['m']/Rv
    eta=4*Kv/num['Q']
    xH=4*Kv/(num['Q']+3*Kv); xHe=4*(4*Kv)/(num['Q']+3*4*Kv)
    return Kv,eta,xH,xHe
print("NS",vals(2.8e30,1.2e4))
print("WD",vals(1.2e30,7.0e6))
RH=num['G']*2.8e30*num['m']/num['Q']; print("R_H",RH,"R_He",4*RH)
print("K in MeV NS", vals(2.8e30,1.2e4)[0]/1.602e-13)
# brute-force check: for WD, scan x, find minimal x where e>=K and e>=4K
import numpy as np
Kv=vals(1.2e30,7.0e6)[0]; xs=np.linspace(0,1,2000001); ev=xs*num['Q']/(4-3*xs)
print("scan xH",xs[np.argmax(ev>=Kv)],"xHe",xs[np.argmax(ev>=4*Kv)])
# check energy conservation of counting: total energy shared = per-atom e * atoms
Nv=1000.;xv=0.3
print(abs(xv*Nv/4*num['Q'] - (Nv*(1-xv)+Nv*xv/4)*xv*num['Q']/(4-3*xv))<1e-20)
