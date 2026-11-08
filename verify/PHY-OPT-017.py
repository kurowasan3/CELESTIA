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
# (4) ratio lambda(theta)/lambda(0) = sqrt(1 - sin^2(theta)/n^2), independent of d and order
th_,ns_,dd_,kk_=sp.symbols('theta n d k',positive=True)
phi_=sp.asin(sp.sin(th_)/ns_)
lam_th=2*ns_*dd_*sp.cos(phi_)/(kk_+sp.Rational(1,2)); lam_0=2*ns_*dd_/(kk_+sp.Rational(1,2))
r_=sp.simplify(lam_th/lam_0); print('ratio =',r_)
assert sp.simplify(r_-sp.sqrt(1-sp.sin(th_)**2/ns_**2))==0
print('check (3): ratio at sin=0.75, n=1.25 ->',float(sp.sqrt(1-sp.Rational(9,16)/n**2)),' 504/630=',504/630)
lim=sp.limit(sp.sqrt(1-sp.sin(th_)**2/ns_**2),th_,sp.pi/2); print('limit theta->90deg:',sp.simplify(lim),' n=1.25 ->',float(lim.subs(ns_,n)),' ->',630*float(lim.subs(ns_,n)),'nm')
assert abs(630*float(lim.subs(ns_,n))-378)<1e-9
# s-polarized exact Airy reflectance: track the k=3 minimum as theta increases (phase rules same as normal incidence for s-pol here)
lams=np.linspace(370e-9,700e-9,200001)
for sth in [0.0,0.5,0.75,0.9,0.99]:
    th=np.arcsin(sth); Rv=R(lams,th,'s'); m_=(Rv[1:-1]<Rv[:-2])&(Rv[1:-1]<=Rv[2:]); mins=lams[1:-1][m_]*1e9
    pred=630*np.sqrt(1-sth**2/nf**2)
    near=mins[np.argmin(abs(mins-pred))]; print('sin th=%.2f predicted %.1f nm, Airy s-pol minimum %.1f nm'%(sth,pred,near))
    assert abs(near-pred)<0.5
# (c) larger n -> smaller shift for the same theta
for nn in [1.25,1.5,2.0,2.5]: print('n=',nn,' ratio at 60deg:',round(np.sqrt(1-np.sin(np.radians(60))**2/nn**2),4))
print('OK')
