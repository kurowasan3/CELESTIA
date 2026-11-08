# O20262-1 verification: constant-voltage Hall sensor and tilted rotation
import sympy as sp, numpy as np
n,q,rho,VE,L,d,h,B,c,v=sp.symbols('n q rho V_E L d h B c v',positive=True)
E=VE/L
csol=sp.solve(sp.Eq(n*q*(q*E/c),E/rho),c)[0]; print('c=',csol)
vd=q*E/csol; VH=sp.simplify(vd*B*d); print('V_H=',VH)
I=sp.simplify(n*q*vd*d*h); print('I=',I, ' check IB/(nqh)=',sp.simplify(I*B/(n*q*h)-VH))
K=float((VH/B).subs({VE:10.0,d:1e-2,L:1e-2,n:5e19,q:1.6e-19,rho:2.0})); print('K=',K)
eps=np.radians(30); Vmax,Vmin=5.0e-6,-27.0e-6
Bh=(Vmax-Vmin)/(2*K*np.cos(eps)); Bv=-(Vmax+Vmin)/(2*K*np.sin(eps))
print('Bh=',Bh,'Bv=',Bv,'dip=',np.degrees(np.arctan2(Bv,Bh)))
# forward model check
th=np.radians(np.arange(0,360,0.01)); thN=np.radians(70)
nvec=np.stack([np.cos(eps)*np.cos(th-thN),np.cos(eps)*np.sin(th-thN),np.sin(eps)*np.ones_like(th)])
Bvec=np.array([Bh,0,-Bv])
V=K*(Bvec@nvec)
print('max',V.max(),'at',np.degrees(th[V.argmax()]),' min',V.min(),'at',np.degrees(th[V.argmin()]))
pos=np.degrees(th[V>0]); print('positive range',pos.min(),pos.max())
print('acos(11/16)=',np.degrees(np.arccos(11/16)))
print('critical eps =90-dip:',90-np.degrees(np.arctan2(Bv,Bh)))
