# O20233-2 rotating observer, fixed source S at distance d from center O
import sympy as sp, numpy as np
r,d,w,th=sp.symbols('r d omega theta',positive=True)
P=sp.Matrix([r*sp.cos(th),r*sp.sin(th)]); S=sp.Matrix([d,0])
v=sp.diff(P,th)*w
SP=P-S
vsp=sp.simplify((v.dot(SP))/sp.sqrt(SP.dot(SP)))
print('v_SP=',vsp)
# distance from O to line SP: |S x P|/|SP|
h=sp.simplify((S[0]*P[1]-S[1]*P[0])/sp.sqrt(SP.dot(SP)))
print('v_SP - omega*h =',sp.simplify(vsp-w*h))
c=sp.symbols('c')
crit=sp.solve(sp.diff((1-c**2)/(r**2+d**2-2*r*d*c),c),c); print('critical cos:',crit)
# numeric scan both cases
for rr,dd in [(5.0,2.0),(5.0,10.0)]:
    t=np.linspace(0,2*np.pi,2000001)
    f=sp.lambdify(th,vsp.subs({r:rr,d:dd,w:2.0}),'numpy'); vv=f(t)
    i=np.argmax(vv); j=np.argmin(vv)
    print(rr,dd,'max vSP',vv[i],'at th',t[i],'min',vv[j],'at',t[j],' expected',2.0*min(rr,dd), np.arccos(min(rr,dd)/max(rr,dd)))
# numeric inverse problem
V=340.0; fmax=700.0; fmin=660.0; wv=2.0; rv=5.0; tau=2*np.pi/3
f0=(fmax+fmin)/2; D=V*(fmax-fmin)/(fmax+fmin)/wv
q=-np.cos(wv*tau/2)
print('f0',f0,'D',D,'q',q,'inside cand d=q r',q*rv,'(needs D=d=',q*rv,') outside d=r/q',rv/q,'(needs D=r=',rv,')')
# simulate actual observed frequency for d=10 and check fmax,fmin,tau
t=np.linspace(0,np.pi*2/wv,2000001); thv=wv*t
f=sp.lambdify(th,vsp.subs({r:rv,d:10.0,w:wv}),'numpy'); vv=f(thv)
fo=f0*(1-vv/V); i=np.argmax(fo); j=np.argmin(fo)
print('sim fmax',fo[i],'fmin',fo[j],'tau(min->max)',t[i]-t[j])
