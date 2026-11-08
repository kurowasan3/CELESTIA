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
# (5) stacked crossed plates (second plate's u along first plate's w)
def stack(lam,d1,d2,dn=0.0090):
    G1=2*np.pi*dn*d1*1000/lam; G2=2*np.pi*dn*d2*1000/lam
    phi=np.pi/4; u=np.array([np.sin(phi),np.cos(phi)]); w=np.array([np.cos(phi),-np.sin(phi)])
    E=np.array([0,1.0],dtype=complex); Eu=np.dot(E,u)*np.exp(-1j*G1); Ew=np.dot(E,w)*np.exp(-1j*G2)
    Eo=Eu*u+Ew*w; return abs(Eo[0])**2
d1=float(D)/0.0090/1000
rng=np.random.default_rng(3)
for _ in range(200):
    lam=rng.uniform(410,700); d2=rng.uniform(0,2*d1)
    assert abs(stack(lam,d1,d2)-np.sin(np.pi*0.0090*(d1-d2)*1000/lam)**2)<1e-12
print('(a) I = I0 sin^2(pi dn (d-d\')/lam) verified')
print('(b) d\'=d: max I over spectrum =',max(stack(l,d1,d1) for l in np.linspace(410,700,301)))
def ndark(X): return [X/k for k in range(1,200) if 410<=X/k<=700]
Xs=np.arange(1,3001,1.0); none=[X for X in Xs if len(ndark(X))==0]
segs=[]; start=none[0]; prev=none[0]
for X in none[1:]:
    if X!=prev+1: segs.append((start,prev)); start=X
    prev=X
segs.append((start,prev)); print('(b) X ranges (nm, integer grid) with no dark line:',segs)
assert segs==[(1.0,409.0),(701.0,819.0)]
for X in [3000,1500,800,300]: print('   X=',X,' dark lines:',[round(v,1) for v in ndark(X)])
# (c) sign symmetry
for lam in [450,520,690]:
    assert abs(stack(lam,d1,d1-50)-stack(lam,d1,d1+50))<1e-12
print('(c) d\'=d-delta and d+delta give identical spectra')
print('OK')
