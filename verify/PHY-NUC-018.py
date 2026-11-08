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
# numbers
dE=2.0; d=0.60e-3
Pperp=2*m*d/R; Pc=np.hypot(dE,Pperp); K=Pc**2/(2*m)*1e3
print('P_perp c=',Pperp,'keV  Pc=',Pc,'keV  K=',K,'eV  alpha=',np.degrees(np.arctan2(Pperp,dE)))
# check K with exact E: electron total energy sqrt(m^2+P^2)-m
print('exact K eV',(np.sqrt(m**2+Pc**2)-m)*1e3)
# symbolic first-order derivation
P,a,M=sp.symbols('P alpha m',positive=True); eps=sp.symbols('epsilon',positive=True)
Et=M+sp.sqrt(M**2+(eps*P)**2)
E1=(Et**2-(eps*P)**2)/(2*(Et-eps*P*sp.cos(a)))
print('E1 series',sp.series(E1,eps,0,2))
