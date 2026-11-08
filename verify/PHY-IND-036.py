# O20132-1 検証：ソレノイドの直流励磁とコンデンサー放電（巻数によらない最大磁場）
from sympy import *
import numpy as np
from scipy.integrate import solve_ivp
mu0,N,A,l,E,r,C,t=symbols('mu0 N A ell E r C t',positive=True)
L=mu0*N**2*A/l
B=lambda I: mu0*N/l*I
# (1)
B1=B(E/r)
# (2) LC: q(0)=CE, I(0)=0 -> I=E*sqrt(C/L) sin(wt)
w=1/sqrt(L*C)
I2=E*sqrt(C/L)
B2=simplify(B(I2)); print('B2=',B2, ' dB2/dN=',simplify(diff(B2,N)))
assert simplify(B2-E*sqrt(mu0*C/(A*l)))==0
# energy check
assert simplify(B2**2/(2*mu0)*A*l - C*E**2/2)==0
t2=simplify(pi/2/w); print('t2=',t2)
# (3)
Ns=solve(Eq(B1,B2),N)[0]; print('N*=',Ns)
assert simplify(Ns - r*sqrt(C*l/(mu0*A)))==0
assert simplify(L.subs(N,Ns)-r**2*C)==0
# (4) ODE solution symbolic
I1=E/r
I=Function('I')
for sgn in [1,-1]:
    sol=dsolve(Eq(I(t).diff(t,2),-I(t)/(L*C)),ics={I(0):I1, I(t).diff(t).subs(t,0):sgn*E/L})
    Ix=sol.rhs
    amp2=simplify((Ix**2+(Ix.diff(t)/w)**2))  # constant
    print('sgn',sgn,'Imax^2=',simplify(amp2))
    assert simplify(amp2-(I1**2+I2**2))==0
# numeric check of times and sign with values
vals={mu0:4e-7*np.pi,N:500,A:1e-3,l:0.2,E:10.0,r:2.0,C:1e-3}
Ln=float(L.subs(vals)); Cn=1e-3; I1n=10/2.0; I2n=float(I2.subs(vals))
for sgn in [1,-1]:
    f=lambda tt,y:[-y[1], y[0]/Cn/Ln*0 + y[1]*0] # placeholder
    # state y=[q,I]: dq/dt=-I, L dI/dt = q/C
    f=lambda tt,y:[-y[1], y[0]/(Cn*Ln)]
    T=2*np.pi*np.sqrt(Ln*Cn)
    s=solve_ivp(f,[0,T],[sgn*Cn*10.0,I1n],max_step=T/20000,rtol=1e-10,atol=1e-12)
    aI=np.abs(s.y[1]); k=int(np.argmax(aI>=aI.max()*(1-1e-9))); tpk=s.t[k]; Ipk=s.y[1][k]
    phi=np.arctan(I2n/I1n); w_=1/np.sqrt(Ln*Cn)
    pred=phi/w_ if sgn>0 else (np.pi-phi)/w_
    print('sgn',sgn,'tpk num',tpk,'pred',pred,'Ipk',Ipk,'pred |I|',np.hypot(I1n,I2n))
    assert abs(tpk-pred)/T<1e-3 and abs(abs(Ipk)-np.hypot(I1n,I2n))<1e-4*abs(Ipk)
    assert (Ipk>0)==(sgn>0)
# (5) new: kappa = Bmax/max(B1,B2) as function of x=N/N*, maximum sqrt2 at N=N*
x=symbols('x',positive=True)
ratio=simplify((B1/B2).subs(N,x*Ns)); print('B1/B2 at N=x N* :',ratio)
assert simplify(ratio-x)==0
k_lo=sqrt(1+x**2); k_hi=sqrt(1+1/x**2)
print('kappa(x<=1)=',k_lo,' kappa(x>=1)=',k_hi,' at x=1:',k_lo.subs(x,1),k_hi.subs(x,1))
xs=np.linspace(0.05,5,1981)  # includes x=1 exactly
kap=np.sqrt(1+xs**2)/np.maximum(1,xs)
print('numeric argmax x=',xs[np.argmax(kap)],' max=',kap.max())
assert abs(xs[np.argmax(kap)]-1)<5e-3 and abs(kap.max()-np.sqrt(2))<1e-4
# phase-quadrature: max |I1 cos + I2 sin| = sqrt(I1^2+I2^2) < I1+I2
tt=np.linspace(0,2*np.pi,100001)
for a,b in [(1,1),(1,0.3),(2,5)]:
    m=np.abs(a*np.cos(tt)+b*np.sin(tt)).max(); print(a,b,'max',m,'hypot',np.hypot(a,b),'sum',a+b)
    assert abs(m-np.hypot(a,b))<1e-6 and m<a+b
print('OK')
