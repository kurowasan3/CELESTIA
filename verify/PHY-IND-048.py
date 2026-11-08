# S0531-1: non-planar loop (cube hexagon) rotating in uniform field; vector area
import numpy as np, sympy as sp
a=1.0
P=np.array([[a,0,0],[a,a,0],[0,a,0],[0,a,a],[0,0,a],[a,0,a]])  # A B C D E F
def vec_area(P):
    S=np.zeros(3)
    for i in range(len(P)):
        S+=0.5*np.cross(P[i],P[(i+1)%len(P)])
    return S
print('vector area of A->B->C->D->E->F',vec_area(P))
# flux by direct surface integral: triangulate fan from centroid & integrate B.n numerically (B uniform -> same); use independent check via Stokes with vector potential A=(1/2)B x r
def flux_line(P,Bv,n=2000):
    tot=0
    for i in range(len(P)):
        p,q=P[i],P[(i+1)%len(P)]
        s=np.linspace(0,1,n+1); mid=(s[:-1]+s[1:])/2
        r=p+np.outer(mid,q-p); dr=(q-p)/n
        Avec=0.5*np.cross(Bv,r)
        tot+=np.sum(Avec@dr)
    return tot
B=1.0; w=1.0
def Rz(th): c,s=np.cos(th),np.sin(th); return np.array([[c,-s,0],[s,c,0],[0,0,1]])
ts=np.linspace(0,2*np.pi,2001)
Phi=np.array([flux_line(P@Rz(t).T,np.array([B,0,0]),200) for t in ts])
print('Phi(0)=',Phi[0],' expect B a^2 =',B*a*a)
emf=-np.gradient(Phi,ts)
print('emf amplitude',emf.max(),' expect sqrt2 w B a^2',np.sqrt(2))
print('|emf| at t=0',abs(emf[0]),' expect wBa^2=1')
print('first t of max |emf|',ts[np.argmax(np.abs(emf[:600]))],' expect pi/4',np.pi/4)
print('Phi formula check',np.max(np.abs(Phi-np.sqrt(2)*B*a*a*np.cos(ts+np.pi/4))))
# inverse problem: random closed non-planar loop
rng=np.random.default_rng(1)
Q=rng.normal(size=(9,3))
S=vec_area(Q)
def amp(Q,axis):
    axis=axis/np.linalg.norm(axis)
    # field perpendicular to axis, rotating loop about axis
    e1=np.cross(axis,[0.3,0.7,0.1]); e1/=np.linalg.norm(e1)
    th=np.linspace(0,2*np.pi,721)
    ph=[]
    for t in th:
        c,s=np.cos(t),np.sin(t); K=np.array([[0,-axis[2],axis[1]],[axis[2],0,-axis[0]],[-axis[1],axis[0],0]])
        Rm=np.eye(3)+s*K+(1-c)*K@K
        ph.append(flux_line(Q@Rm.T,e1*B,40))
    ph=np.array(ph); return np.max(-np.gradient(ph,th))
V=[amp(Q,np.eye(3)[k]) for k in range(3)]
Vmax_formula=np.sqrt(sum(v*v for v in V)/2)
# brute-force max over axes
best=0
for k in range(300):
    ax=rng.normal(size=3); best=max(best,np.linalg.norm(S-np.dot(S,ax/np.linalg.norm(ax))*ax/np.linalg.norm(ax)))
print('V1..V3',V,' formula Vmax',Vmax_formula,' |S|wB',np.linalg.norm(S),' brute(axis sampling)',best)
print('direct amp for axis perp to S:',amp(Q,np.cross(S,[1,0,0])))
print('cube: Vx,Vy,Vz all sqrt2, Vmax sqrt3 ->',np.sqrt((2+2+2)/2))
