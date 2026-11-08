# P52-1: AR-type film, adjacent reflection minima 630,490 nm at normal incidence; minimum shifts to 504 nm at sin(theta)=0.75
import sympy as sp, numpy as np
k=sp.symbols('k',integer=True,nonnegative=True)
# case A (n<nG): 2nd=(k+1/2)lam ; adjacent: (k+1/2)*630=(k+3/2)*490
print('case A k=',sp.solve(sp.Eq((k+sp.Rational(1,2))*630,(k+sp.Rational(3,2))*490),k))
kk=sp.symbols('kk',positive=True)
print('case B k=',sp.solve(sp.Eq(kk*630,(kk+1)*490),kk),'(must be integer)')
two_nd=sp.Rational(7,2)*630; print('2nd =',two_nd,'nm')
cphi=sp.Rational(504,630); sphi=sp.sqrt(1-cphi**2); n=sp.Rational(3,4)/sphi
print('cos phi',cphi,'sin phi',sphi,'n=',n,float(n)); dd=two_nd/(2*n); print('d=',dd,float(dd),'nm')
# transmitted maxima at sin theta=0.75 in 400-700nm: 2nd cos(phi)=(k+1/2)lam
opd=two_nd*cphi; print('2nd cos phi=',opd,[float(opd/(j+sp.Rational(1,2))) for j in range(1,6)])
# exact multi-beam (Airy) reflectance check with Fresnel coefficients (s-pol & p-pol) to confirm minimum location shift
nG=1.5; nf=float(n); d=float(dd)*1e-9
def R(lam,th,pol):
    s0=np.sin(th); c0=np.cos(th); c1=np.sqrt(1-(s0/nf)**2); c2=np.sqrt(1-(s0/nG)**2)
    if pol=='s': r01=(c0-nf*c1)/(c0+nf*c1); r12=(nf*c1-nG*c2)/(nf*c1+nG*c2)
    else: r01=(nf*c0-c1)/(nf*c0+c1); r12=(nG*c1-nf*c2)/(nG*c1+nf*c2)
    b=2*np.pi*nf*d*c1/lam; rr=(r01+r12*np.exp(2j*b))/(1+r01*r12*np.exp(2j*b)); return abs(rr)**2
lams=np.linspace(450e-9,700e-9,250001)
for th in [0,np.arcsin(0.75)]:
    for pol in ['s','p']:
        Rv=R(lams,th,pol); m=(Rv[1:-1]<Rv[:-2])&(Rv[1:-1]<=Rv[2:]); print('theta',round(th,3),pol,'minima nm',np.round(lams[1:-1][m]*1e9,1))
# Snell check: sin theta = n sin phi
print('n sin phi =',float(n*sphi))
