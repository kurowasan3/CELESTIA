# K20213-1 verification: gravity-induced phase in neutron interferometer with lambda/2 contamination
import sympy as sp, numpy as np
from scipy.optimize import minimize_scalar
m,g,h,l,s,lam,al,v0=sp.symbols('m g h ell s lambda_0 alpha v_0',positive=True)
# exact wavelength on DC side
v1=sp.sqrt(v0**2-2*g*s*sp.sin(al)); lam1=h/(m*v1)
dphi_exact=2*sp.pi*l*(1/lam - 1/lam1)   # path ABC: AB at lam0, path ADC: DC at lam1 (AD,BC identical)
dphi_exact=dphi_exact.subs(v0,h/(m*lam))
ser=sp.series(dphi_exact.subs(g,sp.Symbol('eps')*g),sp.Symbol('eps'),0,2).removeO().subs(sp.Symbol('eps'),1)
print('first order dphi =',sp.simplify(ser))
# numbers
M=1.67e-27; G=9.80; H=6.63e-34; A=7.0e-4; L0=2.0e-10
N=M**2*G*A*L0/H**2; q=2*np.pi*N
print('m^2 g/h^2=',M**2*G/H**2,' N=',N,' q=',q)
# exactness check of approximation: v0, g s / v0^2
v=H/(M*L0); print('v0=',v,'(take s=0.02 m) gs/v0^2=',G*0.02/v**2)
# count of alpha in (0,pi/2] with both components constructive: q sin a = 4 pi k
print('both-constructive count:',int(np.floor(q/(4*np.pi))),' lam0-only count:',int(np.floor(q/(2*np.pi))))
# intensity f(phi)=cos^2(phi/2)+cos^2(phi/4) (normalized so f(0)=2)
f=lambda phi: np.cos(phi/2)**2+np.cos(phi/4)**2
res=minimize_scalar(f,bounds=(0.1,2*np.pi),method='bounded')
print('first min phi=',res.x,' 2*arccos(-1/4)=',2*np.arccos(-0.25),' f=',res.fun,' ratio=',res.fun/2,' 7/32=',7/32)
print('sin alpha at first min =',res.x/q,' naive pi/q=',np.pi/q)
# scan to be sure first minimum as alpha increases
a=np.linspace(0,np.pi/2,200001); I=f(q*np.sin(a))
i=np.argmax(np.diff(np.sign(np.diff(I)))>0)+1
print('scan first local min sin a=',np.sin(a[i]),' I/I0=',I[i]/2)
# check q sin a=2pi value: I/I0
print('I/I0 at phi=2pi',f(2*np.pi)/2)
