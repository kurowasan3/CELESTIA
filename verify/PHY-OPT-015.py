# P50-1: three slits, middle slit with phase plate phi (and amplitude k)
import sympy as sp, numpy as np
de,ph,k=sp.symbols('delta phi k',real=True)
amp=sp.exp(sp.I*de)+sp.exp(-sp.I*de)+k*sp.exp(-sp.I*ph)
I=sp.simplify(sp.expand(amp*sp.conjugate(amp)).rewrite(sp.cos))
print('I =',sp.simplify(sp.expand_trig(I)))
I1=sp.simplify(I.subs(k,1))
print('k=1 :',sp.simplify(I1))
for p in [0,sp.pi/2,sp.pi]:
    f=sp.lambdify(de,I1.subs(ph,p)); xs=np.linspace(-2*np.pi,2*np.pi,400001); y=f(xs)
    # local maxima
    m=(y[1:-1]>y[:-2])&(y[1:-1]>=y[2:]); 
    print('phi=',p,' maxima delta/pi:',np.round(xs[1:-1][m]/np.pi,4),' values',np.round(y[1:-1][m],4),' min',round(y.min(),4))
# can center be dark with k=1? |2+e^{-i phi}|>=1
print('min over phi of I(0) k=1:', sp.minimum(I1.subs(de,0),ph,sp.Interval(0,2*sp.pi)))
I2=sp.simplify(I.subs({k:2,ph:sp.pi})); print('k=2,phi=pi:',sp.factor(sp.simplify(I2)))
print('zeros:',sp.solveset(sp.Eq(I2,0),de,sp.Interval(-2*sp.pi,2*sp.pi)))
print('max:', sp.maximum(I2,de,sp.Interval(0,2*sp.pi)))
# thickness: phi = 2pi(n-1)t/lambda
n,t,lam=sp.symbols('n t lambda',positive=True)
print('t for pi/2:',sp.solve(sp.Eq(2*sp.pi*(n-1)*t/lam,sp.pi/2),t),' t for pi:',sp.solve(sp.Eq(2*sp.pi*(n-1)*t/lam,sp.pi),t))
# exact geometry check for phi=pi/2: equal maxima spacing lambda L/(2a)
lv,a,Lv=500e-9,0.2e-3,2.0
x=np.linspace(-2e-2,2e-2,400001)
def r(y): return np.hypot(Lv,x-y)
A=np.exp(2j*np.pi*r(a)/lv)+np.exp(2j*np.pi*r(-a)/lv)+np.exp(2j*np.pi*(r(0)/lv+0.25))
Iv=np.abs(A)**2; m=(Iv[1:-1]>Iv[:-2])&(Iv[1:-1]>=Iv[2:]); pk=x[1:-1][m]
print('exact geom phi=pi/2: peak spacing',np.diff(pk)[:4],'expected',lv*Lv/(2*a),' peak vals',np.round(Iv[1:-1][m][:4],3),' min',round(Iv.min(),3))
# far-field check: a^2/(lambda L) small -> maxima equal
for Lv in [2.0,200.0]:
    x=np.linspace(-2e-2*Lv/2,2e-2*Lv/2,400001)
    def r(y): return np.hypot(Lv,x-y)
    A=np.exp(2j*np.pi*r(a)/lv)+np.exp(2j*np.pi*r(-a)/lv)+np.exp(2j*np.pi*(r(0)/lv+0.25))
    Iv=np.abs(A)**2; m=(Iv[1:-1]>Iv[:-2])&(Iv[1:-1]>=Iv[2:])
    print('L=',Lv,'a^2/(lam L)=',a**2/(lv*Lv),' peak vals',np.round(Iv[1:-1][m][:4],3))
