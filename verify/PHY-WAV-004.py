# M65-1: composite string rho (length l) + 4rho (length l/2), both ends fixed, tension S
import numpy as np, sympy as sp
from scipy.optimize import brentq
from scipy.linalg import eigh
l=1.0; S=1.0; rho=1.0
v1=np.sqrt(S/rho); v2=np.sqrt(S/(4*rho)); l1=l; l2=l/2
print("transit times",l1/v1,l2/v2)
# exact eigen condition from continuity of y and y'
def D(f):
    k1=2*np.pi*f/v1; k2=2*np.pi*f/v2
    return k2*np.sin(k1*l1)*np.cos(k2*l2)+k1*np.cos(k1*l1)*np.sin(k2*l2)
fs=np.linspace(1e-3,2.6,26001); vals=D(fs); roots=[]
for i in range(len(fs)-1):
    if vals[i]*vals[i+1]<0: roots.append(brentq(D,fs[i],fs[i+1]))
print("roots/(v1/4l)",[r/(v1/(4*l)) for r in roots])
# finite-difference independent check
N=1500; L=l1+l2; x=np.linspace(0,L,N+1); dx=x[1]-x[0]
xi=x[1:-1]; mu=np.where(xi<l1,rho,4*rho); mu[np.isclose(xi,l1)]=2.5*rho
K=(np.diag(2*np.ones(N-1))-np.diag(np.ones(N-2),1)-np.diag(np.ones(N-2),-1))*S/dx**2
w2,vec=eigh(K,np.diag(mu))
f_fd=np.sqrt(w2[:6])/(2*np.pi); print("FD f/(v1/4l)",f_fd/(v1/(4*l)))
for j in range(4):
    y=vec[:,j]; iB=np.argmin(abs(xi-l1))
    A1=np.max(abs(y[xi<l1])); A2=np.max(abs(y[xi>l1]))
    print("mode",j+1,"yB/A1=%.3f"%(y[iB]/A1 if A1 else 0),"A2/A1=%.3f"%(A2/A1))
# symbolic: equal transit -> tan(theta)(1/k1+1/k2)=0
th=sp.symbols('theta')
print(sp.simplify((sp.Rational(1,1)*sp.tan(th)+sp.tan(th)/2)))
# (5) M_max = 16 rho l^2 f0^2 / g
