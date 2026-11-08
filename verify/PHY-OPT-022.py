# K20193-1: multiple-beam interference in a film (Fabry-Perot)
import numpy as np, sympy as sp
def trans(p1,p2,delta,K=4000):
    # amplitude: q1 q2' * sum (p1 p2 e^{i delta})^k ; with q1q1'=1-p1^2 ... use t=sqrt products
    t=np.sqrt((1-p1**2)*(1-p2**2))
    s=sum((p1*p2*np.exp(1j*delta))**kk for kk in range(K))
    return abs(t*s)**2
def refl(p1,p2,delta,K=4000):
    # first reflection at A from outside: -p1 ; others: (1-p1^2) p2 e^{i delta} sum (p1p2 e^{i delta})^k
    s=sum((p1*p2*np.exp(1j*delta))**kk for kk in range(K))
    return abs(-p1+(1-p1**2)*p2*np.exp(1j*delta)*s)**2
p=np.sqrt(0.8)
for dl in [0,0.3,np.pi]:
    formula=(1-p**2)**2/(1-2*p**2*np.cos(dl)+p**4)
    print('delta',dl,'series',trans(p,p,dl),'formula',formula,'R+T',trans(p,p,dl)+refl(p,p,dl))
print('min/max',trans(p,p,np.pi)/trans(p,p,0),'=1/81?',1/81)
# data inversion
l1,l2=480e-9,600e-9
twonD=1/(1/l1-1/l2); print('2nD nm',twonD*1e9,'orders',twonD/l1,twonD/l2,'D(n=1.5) nm',twonD/3*1e9)
# FWHM near 600 nm exact
lam=np.linspace(580e-9,620e-9,400001)
I=np.array([(1-p**2)**2/(1-2*p**2*np.cos(2*np.pi*twonD/l)+p**4) for l in lam])
sel=lam[I>=0.5]; print('FWHM nm',(sel.max()-sel.min())*1e9)
dl=(1-p**2)/p; print('approx FWHM',dl*(600e-9)**2/(2*np.pi*twonD)*2*1e9)
# reflection at resonance
print('R at resonance',refl(p,p,0))
# unequal
p1,p2=np.sqrt(0.8),np.sqrt(0.2)
print('Tpeak unequal',trans(p1,p2,0),'formula',(1-p1**2)*(1-p2**2)/(1-p1*p2)**2,'R',refl(p1,p2,0))
x,y=sp.symbols('x y',positive=True)
print(sp.factor(1-(1-x**2)*(1-y**2)/(1-x*y)**2))
