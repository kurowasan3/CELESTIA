# O20263-1 verification: exact ray tracing through a glass hemisphere (flat face toward light)
import sympy as sp, numpy as np
R,n,y=sp.symbols('R n y',positive=True)
th=sp.asin(y/R); ph=sp.asin(n*y/R)
CF=R*sp.sin(ph)/sp.sin(ph-th)
print('paraxial:',sp.limit(CF,y,0))
CFc=sp.limit(CF,y,R/n,'-'); print('at y->R/n:',sp.simplify(CFc))
v=CF.subs({n:sp.sqrt(3),y:R/2}); print('n=sqrt3,y=R/2:',sp.nsimplify(sp.simplify(v)))
print('paraxial n=sqrt3:',sp.radsimp(sp.limit(CF,y,0).subs(n,sp.sqrt(3))))
# numerical ray tracing independent: vector Snell
def trace(yv,nv,nout=1.0,Rv=1.0):
    xA=np.sqrt(Rv**2-yv**2); P=np.array([xA,yv]); nrm=P/Rv; d=np.array([1.0,0])
    mu=nv/nout; cosi=d@nrm; s2=mu**2*(1-cosi**2)
    if s2>1: return None
    t=mu*d+(np.sqrt(1-s2)-mu*cosi)*nrm
    s=-P[1]/t[1]; return P[0]+s*t[0], t
for yv in [1e-4,0.3,0.5,0.57]:
    print(yv, trace(yv,np.sqrt(3))[0], float(CF.subs({R:1,n:sp.sqrt(3),y:yv})))
print('TIR at 0.58:',trace(0.58,np.sqrt(3)))
# monotonic decrease check
ys=np.linspace(1e-4,1/np.sqrt(3)-1e-6,2000); xs=[trace(t,np.sqrt(3))[0] for t in ys]
print('monotone decreasing:',np.all(np.diff(xs)<0))
# blur at paraxial focus screen for y=R/2
xp=(3+np.sqrt(3))/2; x1,t=trace(0.5,np.sqrt(3)); yscreen=0+(xp-x1)*t[1]/t[0]
print('screen height',yscreen,' expected -(sqrt3-1)/2=',-(np.sqrt(3)-1)/2)
# liquid nA=1.2 (n=sqrt3): paraxial CF and yc
nA=sp.symbols('n_A',positive=True)
CFl=sp.limit(CF.subs(n,n/nA),y,0); print('liquid paraxial CF:',sp.simplify(CFl),' VF:',sp.simplify(CFl-R))
print('numeric liquid check',trace(1e-5,np.sqrt(3),1.2)[0], float((sp.sqrt(3)/(sp.sqrt(3)-1.2))))
