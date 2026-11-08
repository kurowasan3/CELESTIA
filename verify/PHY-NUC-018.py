import sympy as sp, numpy as np
# units: c=1, energies in keV, m=511
m=511.0
def exact(P,al):
    Et=m+np.sqrt(m**2+P**2)   # electron moving with momentum P, positron at rest
    E1=(Et**2-P**2)/(2*(Et-P*np.cos(al)))
    p1=E1*np.array([np.cos(al),np.sin(al)])
    p2=np.array([P,0])-p1; E2=np.linalg.norm(p2)
    assert abs(E1+E2-Et)<1e-9
    # deviation from back-to-back
    cosd=np.dot(-p1,p2)/(np.linalg.norm(p1)*E2); delta=np.arccos(np.clip(cosd,-1,1))
    return E1,E2,delta
for P,al in [(2.5,0.3),(2.5,1.2),(5,2.0),(1,np.pi/2)]:
    E1,E2,d=exact(P,al)
    print(f'P={P} al={al:.2f}: E1-E2={E1-E2:.5f} approx {P*np.cos(al):.5f}; E1={E1:.4f} approx {m+P*np.cos(al)/2:.4f}; delta={d:.6e} approx {P*np.sin(al)/m:.6e}')
# chord distance: two points on circle radius R, central angle pi-delta -> distance from center R sin(delta/2)
R=0.40
th=np.linspace(0,1,5)
for dl in [1e-3,5e-3]:
    A=R*np.array([1,0]); B=R*np.array([np.cos(np.pi-dl),np.sin(np.pi-dl)])
    # distance from origin to line AB
    u_=B-A; dist=abs(u_[0]*(-A[1])-u_[1]*(-A[0]))/np.linalg.norm(u_)
    print('chord dist',dist,' R*delta/2',R*dl/2)
# symbolic first-order derivation
P,a,M=sp.symbols('P alpha m',positive=True); eps=sp.symbols('epsilon',positive=True)
Et=M+sp.sqrt(M**2+(eps*P)**2)
E1=(Et**2-(eps*P)**2)/(2*(Et-eps*P*sp.cos(a)))
print('E1 series',sp.series(E1,eps,0,2))

# ---- new (4): inversion formulas, d_max(K), d proportional to R ----
# exact kinematics, positron at rest, electron momentum P at angle alpha (relative to photon 1)
def measure(P,al,R):
    E1,E2,dl=exact(P,al)
    A=R*np.array([1,0]); B=R*np.array([np.cos(np.pi-dl),np.sin(np.pi-dl)])
    u_=B-A; dist=abs(u_[0]*(-A[1])-u_[1]*(-A[0]))/np.linalg.norm(u_)
    return E1-E2,dist
for P_,al_ in [(2.5,0.66),(4.0,1.9)]:
    dE_,dd=measure(P_,al_,0.40)
    Prec=np.hypot(dE_,2*m*dd/0.40); print(' true P,cos',P_,np.cos(al_),' reconstructed',Prec,dE_/Prec)
Kev=6.0e-3  # keV
Pk=np.sqrt(2*m*Kev)   # keV (c=1)
for R_ in [0.2,0.4,0.8]:
    dmax=max(measure(Pk,a_,R_)[1] for a_ in np.linspace(0,np.pi,2001))
    print(' R=',R_,' dmax num [mm]',dmax*1e3,' formula R sqrt(K/2mc2) [mm]',R_*np.sqrt(Kev/(2*m))*1e3)
Ks,Rs,mc2s=sp.symbols('K R mc2',positive=True)
print('d_max symbolic check:',sp.simplify(Rs*sp.sqrt(2*Ks*mc2s)/(2*mc2s)-Rs*sp.sqrt(Ks/(2*mc2s))))
