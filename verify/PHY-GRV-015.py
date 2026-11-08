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

# ---- (4) frame co-rotating with V: E moves on circle b with speed b*omega; component along VE = b w sin(SEV) <= a w
ang=np.linspace(0,2*np.pi,200001)
A=1.0*aa; B=1.0; W=1.0
Vp=np.array([A,0.0])                       # V fixed
Ex=B*np.cos(ang); Ey=B*np.sin(ang)          # E positions
vx=B*W*np.sin(ang); vy=-B*W*np.cos(ang)     # E velocity (angular velocity -W)
lx,ly=Ex-Vp[0],Ey-Vp[1]; L=np.hypot(lx,ly)
comp=(vx*lx+vy*ly)/L                        # d(dist)/dt
# sin(SEV): angle at E between ES and EV
cosSEV=((-Ex)*(-lx)+(-Ey)*(-ly))/(B*L); sinSEV=np.sqrt(1-cosSEV**2)
assert np.allclose(np.abs(comp),B*W*sinSEV,atol=1e-9)
print('(4) |comp| = b w sin(SEV) OK; max |comp| =',np.abs(comp).max(),' a*w =',A*W,' b*w =',B*W)
assert abs(np.abs(comp).max()-A*W)<1e-6
# compare with analytic rate formula from (1)
rate_num=sp.lambdify(phi,rate.subs({a:A,b:B,w:W}))
assert np.allclose(np.abs(comp),np.abs(rate_num(ang)),atol=1e-9)
print('(4) agrees with (1) formula')
