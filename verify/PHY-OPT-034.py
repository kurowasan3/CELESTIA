# O20253-1 verification: wedge / parallel films with phase-shift arrangements
import sympy as sp, numpy as np
lam0,L,D,n,F0=sp.symbols('lambda0 L D n F0',positive=True)
# (1) air: one pi shift -> dark condition 2t = m lam ; t = D x/L ; spacing
dx_air = lam0*L/(2*D)
# (2) liquid n (n1<n<n2): both pi -> bright 2 n t = m lam, spacing
Dp=sp.symbols('Dp',positive=True)
dx_liq = lam0*L/(2*n*Dp)
sol=sp.solve(sp.Eq(dx_liq,dx_air),Dp)[0]; print('D\'=',sol)
k=sp.simplify(F0/(D-sol)); print('k=',k)
# II numeric scan
lmin,lmax=400e-9,760e-9
def reinforced(d,mode,nl=1.60):
    out=[]
    for m in range(0,40):
        if mode=='liq':
            if m==0: continue
            lam=2*nl*d/m
        else:
            lam=4*d/(2*m+1)
        if lmin<=lam<=lmax: out.append(lam)
    return out
print('F=0 liquid:',[round(x*1e9,1) for x in reinforced(400e-9,'liq')])
def windows(mode):
    ds=np.linspace(1e-10,400e-9,400001); empty=[len(reinforced(d,mode))==0 for d in ds]
    res=[];start=None
    for d,e in zip(ds,empty):
        if e and start is None: start=d
        if (not e) and start is not None: res.append((start,d)); start=None
    if start is not None: res.append((start,ds[-1]))
    return [(round(a*1e9,2),round(b*1e9,2)) for a,b in res]
print('liquid windows nm:',windows('liq'))
print('air windows nm:',windows('air'))
# general ratio condition: liquid consecutive ratio (m+1)/m max 2 ; air (2m+3)/(2m+1) max 3
r=sp.symbols('r',positive=True)
print('liquid window between m=1,2 exists iff 2>r ; air iff 3>r ; r here=',lmax/lmin)
