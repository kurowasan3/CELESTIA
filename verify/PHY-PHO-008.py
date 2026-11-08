# T20053-1 falling-atom double slit: fringe spacing = h*T/(m d)
import sympy as sp, numpy as np
H,l,g,hP,m,d,x,v1,v2=sp.symbols('H l g h m d x v1 v2',positive=True)
V1=sp.sqrt(2*g*H); V2=sp.sqrt(2*g*(H+l))
T=sp.simplify((V2-V1)/g)
print('T*(v1+v2)/2 == l :',sp.simplify(T*(V1+V2)/2-l)==0)
dx=hP*T/(m*d); dx2=2*hP*l/(m*d*(V1+V2))
print('two forms equal:',sp.simplify(dx-dx2)==0)
ratio=sp.simplify(dx/(hP/(m*V1)*l/d)); print('ratio to naive:',ratio)
print('limit l->0:',sp.limit(ratio,l,0))
# numeric (4)
vals={H:0.050,l:0.85,g:9.8,hP:6.6e-34,m:3.3e-26,d:6.0e-6}
print('Delta x =',float(dx.subs(vals)),'m')
# independent check by explicit phase computation (semi-classical): phase = action S/hbar along classical paths
# transverse: free motion from slit at x=+-d/2 to screen point X in time T -> phase m (X-xs)^2/(2 hbar T)
hb=6.6e-34/(2*np.pi); mm=3.3e-26; dd=6e-6
Tn=float(T.subs(vals))
X=np.linspace(-5e-3,5e-3,200001)
ph=mm*((X-dd/2)**2-(X+dd/2)**2)/(2*hb*Tn)
I=np.cos(ph/2)**2
pk=X[1:-1][(I[1:-1]>I[:-2])&(I[1:-1]>I[2:])]
print('numeric fringe spacing',np.mean(np.diff(pk)))
assert abs(np.mean(np.diff(pk))/float(dx.subs(vals))-1)<1e-3
# horizontal version: vertical slits separated by d, horizontal speed v, distance L
v,L=sp.symbols('v L',positive=True)
TT=L/v
# arrival z = uz*TT - g TT^2/2 ; phase diff per unit uz -> spacing in z equals h TT/(m d), shifted by -g TT^2/2
uz=sp.symbols('u_z'); z=sp.symbols('z')
uzsol=sp.solve(sp.Eq(z,uz*TT-g*TT**2/2),uz)[0]
# phase difference 2*pi*d*m*uz/h ; spacing in z where it changes by 2pi
spacing=sp.simplify(sp.solve(sp.Eq(2*sp.pi*d*m*sp.diff(uzsol,z)*z/hP,2*sp.pi),z)[0]); print('horizontal spacing:',spacing)
print('shift (center uz=0):',sp.simplify(sp.solve(sp.Eq(uzsol,0),z)[0]))
print('OK')
