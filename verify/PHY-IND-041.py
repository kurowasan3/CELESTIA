# Q176-1 rotating ring in Earth's horizontal field, compass at centre, with self-inductance
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
a,R,L,w,Bh,mu,t=sp.symbols('a R L omega B_h mu0 t',positive=True)
phi=sp.atan(w*L/R); Z=sp.sqrt(R**2+w**2*L**2); S=sp.pi*a**2
I=Bh*S*w/Z*sp.sin(w*t-phi)
# check ODE L I' + R I = Bh S w sin(wt)
print('ode check',sp.simplify(sp.expand_trig(L*sp.diff(I,t)+R*I-Bh*S*w*sp.sin(w*t))))
T=2*sp.pi/w
Bx=sp.simplify(w/(2*sp.pi)*sp.integrate(mu*I/(2*a)*sp.cos(w*t),(t,0,T)))
By=sp.simplify(w/(2*sp.pi)*sp.integrate(mu*I/(2*a)*sp.sin(w*t),(t,0,T)))
print('Bx=',sp.simplify(Bx),' By=',sp.simplify(By))
k=sp.symbols('k',positive=True)  # k=mu0 pi a/(4L)
x=sp.symbols('x',positive=True)  # x=wL/R
Bx_k=sp.simplify((Bx/Bh).subs(mu,4*L*k/(sp.pi*a)).subs(w,x*R/L)); By_k=sp.simplify((By/Bh).subs(mu,4*L*k/(sp.pi*a)).subs(w,x*R/L))
print('Bx/Bh=',sp.factor(Bx_k),' By/Bh=',sp.factor(By_k))
tan=sp.simplify(By_k/(1+Bx_k)); print('tan theta=',sp.factor(tan))
xs=sp.solve(sp.diff(tan,x),x); print('x*=',xs,' tanmax=',[sp.simplify(tan.subs(x,s)) for s in xs])
print('limit x->oo: Bx/Bh',sp.limit(Bx_k,x,sp.oo),' By',sp.limit(By_k,x,sp.oo))
print('L->0 tan:',sp.simplify(sp.limit((By/(Bh+Bx)),L,0)))
# numeric ODE check with values
av,Rv,Lv,wv,Bhv,muv=0.10,2.0e-4,6.0e-7,2*np.pi*50,3.0e-5,4e-7*np.pi
f=lambda tt,y:[(Bhv*np.pi*av**2*wv*np.sin(wv*tt)-Rv*y[0])/Lv]
Tp=2*np.pi/wv; so=solve_ivp(f,[0,60*Tp],[0.0],max_step=Tp/400,rtol=1e-10,atol=1e-14,dense_output=True)
tt=np.linspace(50*Tp,60*Tp,200001); Ii=so.sol(tt)[0]
bx=np.mean(muv*Ii/(2*av)*np.cos(wv*tt)); by=np.mean(muv*Ii/(2*av)*np.sin(wv*tt))
kv=muv*np.pi*av/(4*Lv); xv=wv*Lv/Rv
print('num bx/Bh',bx/Bhv,' formula',-kv*xv**2/(1+xv**2),' by/Bh',by/Bhv,' formula',kv*xv/(1+xv**2))
# P176-original-like numbers for sanity of (1): tan2deg
print('R for 2deg (orig-like)', muv*np.pi*0.1*(2*np.pi*10)/(4*np.tan(np.radians(2))))
