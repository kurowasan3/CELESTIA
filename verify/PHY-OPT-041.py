# S0411-1 Young double slit with polarizers on each slit, elliptically polarized light
import numpy as np, sympy as sp
def pattern(a,b,theta,phis):
    t=np.linspace(0,2*np.pi,4000,endpoint=False)  # psi = omega t
    out=[]
    eth=np.array([np.cos(theta),np.sin(theta)])
    for phi in phis:
        # slit1 polarizer along x ; slit2 along e_theta, its light delayed by phase phi
        E1x=a*np.cos(t); E1=np.array([E1x,0*t])
        psi2=t-phi
        Eth=a*np.cos(theta)*np.cos(psi2)+b*np.sin(theta)*np.sin(psi2)
        E2=np.outer(eth,Eth)
        E=E1+E2
        out.append(np.mean(E[0]**2+E[1]**2)*2)  # x2: amplitude^2 normalization
    return np.array(out)
phis=np.linspace(-np.pi,np.pi,3601)
def analyse(a,b,theta):
    I=pattern(a,b,theta,phis); k=np.argmax(I)
    V=(I.max()-I.min())/(I.max()+I.min())
    return phis[k],V
# circular light
for th in [np.pi/6,np.pi/4,np.pi/3]:
    phimax,V=analyse(1,1,th)
    print('circular theta=%.3f: bright at phi=%.4f (theory %.4f), V=%.4f (theory %.4f)'%(th,phimax,-th,V,np.cos(th)))
# elliptic b=sqrt3 a, theta=45
phimax,V=analyse(1,np.sqrt(3),np.pi/4)
print('b/a=sqrt3, theta=45: phi_bright=%.4f (theory -pi/3=%.4f), V=%.4f (theory 2/3)'%(phimax,-np.pi/3,V))
# without polarizer on slit2 check intensity sums with a,b (original-type): skip
# symbolic formula
a,b,th,phi,psi=sp.symbols('a b theta phi psi',real=True)
E2=a*sp.cos(th)*sp.cos(psi-phi)+b*sp.sin(th)*sp.sin(psi-phi)
Ex=a*sp.cos(psi)+E2*sp.cos(th); Ey=E2*sp.sin(th)
I=sp.integrate(sp.expand(Ex**2+Ey**2),(psi,0,2*sp.pi))/sp.pi
I=sp.simplify(sp.expand_trig(I))
A2=a**2*sp.cos(th)**2+b**2*sp.sin(th)**2
beta=sp.atan(b/a*sp.tan(th))
test=a**2+A2+2*a*sp.sqrt(A2)*sp.cos(th)*sp.cos(phi+beta)
import random
for _ in range(5):
    vals={a:random.uniform(0.5,2),b:random.uniform(0.5,2),th:random.uniform(-1.4,1.4),phi:random.uniform(-3,3)}
    print('symbolic check',float(I.subs(vals)),float(test.subs(vals)))
