# N61-1 spectroscopic binary with systemic velocity, unequal masses
import numpy as np, sympy as sp
G=6.67e-11; c=3.0e8; Ry=1.10e7
lam0=16/(3*Ry)   # H-beta: 1/lam = R(1/4-1/16)
print('lam0 =',lam0)
# choose truth
M1,M2=3.0e30,1.5e30; a=2.0e10; V=2.0e4
T=2*np.pi*np.sqrt(a**3/(G*(M1+M2))); w=2*np.pi/T
a1=a*M2/(M1+M2); a2=a*M1/(M1+M2); v1=w*a1; v2=w*a2
t=np.linspace(0,T,400001)
# star1 at angle wt, star2 opposite; observer along +x far away (receding = +x). radial vel = x-velocity
vr1=V+(-v1*np.sin(w*t)); vr2=V+(v2*np.sin(w*t))
l1=lam0*(1+vr1/c); l2=lam0*(1+vr2/c)   # moving source, observer at rest
sep=np.abs(l1-l2)
merge=t[np.where(np.diff(np.sign(l1-l2))!=0)[0]]; print('merge times',merge[:3],'T',T)
tau=merge[1]-merge[0] if merge[0]>1 else merge[1]
lamA=lam0*(1+V/c); i=np.argmax(sep); lB=min(l1[i],l2[i]); lC=max(l1[i],l2[i])
# inversion formulas
Tm=2*tau; Vm=c*(lamA-lam0)/lam0; vB=c*(lamA-lB)/lam0; vC=c*(lC-lamA)/lam0
MB=vC*(vB+vC)**2*Tm/(2*np.pi*G); MC=vB*(vB+vC)**2*Tm/(2*np.pi*G)
print('T',Tm,T,' V',Vm,V,' vB vC',vB,vC,'(true',v1,v2,')')
print('MB MC',MB,MC,' true',M1,M2, ' separation', (vB+vC)*Tm/(2*np.pi), a)
# symbolic check of mass formula
vb,vc,TT,GG=sp.symbols('v_B v_C T G',positive=True)
aa=(vb+vc)*TT/(2*sp.pi); aB=vb*TT/(2*sp.pi)
MCs=vb**2*aa**2/(GG*aB)   # M_B v_B^2/a_B = G M_B M_C/a^2
MBs=vc**2*aa**2/(GG*(vc*TT/(2*sp.pi)))
print(sp.simplify(MBs), sp.simplify(MCs), sp.simplify(MBs+MCs - (vb+vc)**3*TT/(2*sp.pi*GG)))
# Kepler 3rd check
print('Kepler', sp.simplify(4*sp.pi**2*aa**3/(GG*TT**2) - (MBs+MCs)))
