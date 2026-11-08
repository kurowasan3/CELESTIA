# Q117-1: uniformly charged rod; ellipses with foci at the ends are equipotentials; bead on elliptic wire
import numpy as np, sympy as sp
from scipy.integrate import quad
k=lam=1.0; c=1.0; a=1.6; b=np.sqrt(a**2-c**2)
def E(P):
    fx=quad(lambda x:k*lam*(P[0]-x)/((P[0]-x)**2+P[1]**2)**1.5,-c,c,epsabs=1e-12,limit=200)[0]
    fy=quad(lambda x:k*lam*P[1]/((P[0]-x)**2+P[1]**2)**1.5,-c,c,epsabs=1e-12,limit=200)[0]
    return np.array([fx,fy])
def V(P): return quad(lambda x:k*lam/np.hypot(P[0]-x,P[1]),-c,c,epsabs=1e-12,limit=200)[0]
# bisector & arc formula at arbitrary C
C=np.array([0.5,0.8]); A=np.array([-c,0]); B=np.array([c,0])
ua=(A-C)/np.linalg.norm(A-C); ub=(B-C)/np.linalg.norm(B-C); bis=-(ua+ub); bis/=np.linalg.norm(bis)
EC=E(C); gam=np.arccos(np.dot(ua,ub))
print('E dir',EC/np.linalg.norm(EC),' bisector(away)',bis,' |E|',np.linalg.norm(EC),' 2k lam sin(gam/2)/d',2*k*lam*np.sin(gam/2)/C[1])
# equipotential on ellipse
ts=np.linspace(0.01,np.pi-0.01,7)
print('V on ellipse',[round(V(np.array([a*np.cos(t),b*np.sin(t)])),10) for t in ts],' theory',k*lam*np.log((a+c)/(a-c)))
# E at minor/major ends
print('E minor',np.linalg.norm(E(np.array([0,b]))),' theory',2*k*lam*c/(a*b))
print('E major',np.linalg.norm(E(np.array([a,0]))),' theory',2*k*lam*c/b**2)
# bead: charge -q, mass m; speed constant; normal force
m=1.0;q=1.0
v2=2*k*lam*q*c*a/(m*b**2)    # wire force zero at minor end
rc_minor=a**2/b; rc_major=b**2/a
print('force minor',m*v2/rc_minor-q*np.linalg.norm(E(np.array([0,b]))))
Nmaj=m*v2/rc_major-q*np.linalg.norm(E(np.array([a,0]))); print('N major',Nmaj,' theory',2*k*lam*q*c**3/b**4)
# simulate bead constrained on ellipse to confirm constant speed (tangential force zero)
def tang(t):
    P=np.array([a*np.cos(t),b*np.sin(t)]); T=np.array([-a*np.sin(t),b*np.cos(t)]); T/=np.linalg.norm(T)
    return np.dot(-q*E(P),T)
print('max tangential force',max(abs(tang(t)) for t in np.linspace(0,2*np.pi,25)))
