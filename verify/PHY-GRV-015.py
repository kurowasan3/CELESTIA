import sympy as sp, numpy as np
a,b,w,phi=sp.symbols('a b omega phi',positive=True)
d=sp.sqrt(a**2+b**2-2*a*b*sp.cos(phi))
rate=a*b*w*sp.sin(phi)/d   # w = omega_V - omega_E
c=sp.symbols('c')
num=sp.numer(sp.together(sp.diff(rate,phi)))
print(sp.factor(sp.simplify(num.subs(sp.cos(phi),c).subs(sp.sin(phi)**2,1-c**2))))
print('max', sp.simplify(rate.subs(phi,sp.acos(a/b))))
# numeric brute force with Kepler
G_M=1.0; bb=1.0; aa=0.723
wV=np.sqrt(1/aa**3); wE=1.0
t=np.linspace(0,50,2000001)
xV=aa*np.cos(wV*t); yV=aa*np.sin(wV*t); xE=np.cos(t); yE=np.sin(t)
dist=np.hypot(xV-xE,yV-yE); dd=np.gradient(dist,t)
vE=29.8
print('max rate km/s',np.abs(dd).max()*vE, 'formula', vE*(1/np.sqrt(aa)-aa))
i=np.argmax(dd); 
# elongation angle at max
ex,ey=xV[i]-xE[i],yV[i]-yE[i]; sx,sy=-xE[i],-yE[i]
print('elong deg',np.degrees(np.arccos((ex*sx+ey*sy)/np.hypot(ex,ey)/np.hypot(sx,sy))),np.degrees(np.arcsin(aa)))
print('doppler kHz',2*vE*1e3*(1/np.sqrt(aa)-aa)/3.0e8*1.0e9/1e3)
