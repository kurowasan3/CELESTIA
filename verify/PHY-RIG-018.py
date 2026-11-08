# Q069-1: rod / plate hung from one point by strings; tension proportional to string length
import numpy as np, sympy as sp
def solve_rod(L,h1,h2,ag,W=1.0):
    # A=(0,0), B=(L,0) rotated; find P position from string lengths then rotate so G below P
    # P in rod frame: |PA|=h1,|PB|=h2
    px=(h1**2-h2**2+L**2)/(2*L); py=np.sqrt(h1**2-px**2)
    A=np.array([0,0.]);B=np.array([L,0.]);P=np.array([px,py]);G=np.array([ag,0.])
    u=(P-G)/np.linalg.norm(P-G)                 # this direction must be 'up'
    # rotate so that u -> (0,1)
    c,s=u[1],u[0]; R=np.array([[c,-s],[s,c]])
    A,B,P,G=[R@v for v in (A,B,P,G)]
    assert abs((P-G)[0])<1e-12
    # unknown tensions: T1*(P-A)/h1 + T2*(P-B)/h2 + (0,-W)=0
    Mx=np.column_stack([(P-A)/h1,(P-B)/h2]); T=np.linalg.solve(Mx,[0,W])
    return T,np.linalg.norm(P-G)
T,PG=solve_rod(0.6,0.5,0.4,0.3); print('uniform: T',T,' T1/T2',T[0]/T[1],' h1/h2',0.5/0.4,' formula',[h/(2*PG) for h in (0.5,0.4)],' PG',PG,np.sqrt((0.25+0.16)/2-0.09))
T,PG=solve_rod(0.6,0.5,0.4,0.2); print('AG:GB=1:2: T',T,' T1/T2',T[0]/T[1],' 2h1/h2',2*0.5/0.4,' PG',PG, np.sqrt((0.4*0.25+0.2*0.16)/0.6-0.08))
# 3D equilateral plate side a, strings h1,h2,h3
def plate(a,h):
    V=np.array([[0,0,0],[a,0,0],[a/2,a*np.sqrt(3)/2,0]])
    # find P with |P-Vi|=hi (z>0)
    from scipy.optimize import fsolve
    P=fsolve(lambda p:[np.linalg.norm(p-V[i])-h[i] for i in range(3)],[a/2,a/3,1.0])
    G=V.mean(axis=0); up=(P-G)/np.linalg.norm(P-G)
    Mx=np.column_stack([(P-V[i])/h[i] for i in range(3)]); T=np.linalg.solve(Mx,up*1.0)  # W=1 acting -up
    return T,np.linalg.norm(P-G)
h=np.array([0.5,0.6,0.7]); a=0.4
T,PG=plate(a,h); print('plate T',T,' T/h',T/h,' formula PG',PG,np.sqrt((h@h-a*a)/3),' T formula',h/(3*PG))
