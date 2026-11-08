# K20233-2: birefringent plate between polarizers (channeled spectrum), stacked crossed plates
import numpy as np, sympy as sp
def out(phi,Gam,parallel=False,I0=1.0):
    # input along y; plate axes: u (slow, angle phi from y), w (perp)
    u=np.array([np.sin(phi),np.cos(phi)]); w=np.array([np.cos(phi),-np.sin(phi)])
    E=np.array([0,1.0],dtype=complex)
    Eu=np.dot(E,u)*np.exp(-1j*Gam); Ew=np.dot(E,w)
    Eo=Eu*u+Ew*w
    ax=np.array([0,1.0]) if parallel else np.array([1.0,0])
    return abs(np.dot(Eo,ax))**2*I0
for phi in [0.3,np.pi/4,np.pi/8]:
    for G in [0.7,2.0,np.pi]:
        assert abs(out(phi,G)-np.sin(2*phi)**2*np.sin(G/2)**2)<1e-12
        assert abs(out(phi,G,True)-(1-np.sin(2*phi)**2*np.sin(G/2)**2))<1e-12
print('formula ok')
# (2) adjacent darks 600,500 -> D=Dn*d
m=sp.symbols('m',positive=True,integer=True)
msol=sp.solve(sp.Eq(m*600,(m+1)*500),m); print('m=',msol)
D=msol[0]*600; print('Dn d =',D,'nm ; d=',D/0.0090/1000,'um')
print('darks in [410,700]:',[D/k for k in range(1,20) if 410<=D/k<=700])
print('brights:',[D/(k+0.5) for k in range(0,20) if 410<=D/(k+0.5)<=700])
# (4) stacked crossed plates: effective Dn*|d1-d2|=X; exactly one dark at 500 and no full-bright in range
for X in [500,1000,1500,2000]:
    darks=[X/k for k in range(1,40) if 410<=X/k<=700]; brights=[X/(k+0.5) for k in range(0,40) if 410<=X/(k+0.5)<=700]
    print(X,'darks',darks,'brights',brights)
d1=float(D)/0.0090/1000; dd=500/0.0090/1000
print('d2 options',d1-dd,d1+dd)
# numeric check with stacked plates: plate1 slow axis at pi/4, plate2 slow axis at -pi/4 (crossed)
def stack(lam,d1,d2,dn=0.0090):
    G1=2*np.pi*dn*d1*1000/lam; G2=2*np.pi*dn*d2*1000/lam
    # Jones in (u,w) basis of plate1: plate2 swaps slow axis -> retardation -G2
    phi=np.pi/4; u=np.array([np.sin(phi),np.cos(phi)]); w=np.array([np.cos(phi),-np.sin(phi)])
    E=np.array([0,1.0],dtype=complex); Eu=np.dot(E,u)*np.exp(-1j*G1); Ew=np.dot(E,w)*np.exp(-1j*G2)
    Eo=Eu*u+Ew*w; return abs(Eo[0])**2
for d2 in [d1-dd,d1+dd]:
    lams=np.linspace(410,700,29001); I=np.array([stack(l,d1,d2) for l in lams])
    print('d2=%.1f'%d2,' min I at',lams[I.argmin()],I.min(),' max I',I.max())
