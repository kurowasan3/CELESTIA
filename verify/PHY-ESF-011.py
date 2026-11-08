# J102-1: deflection plates with sinusoidal voltage, transit-time effect
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
e,m,d,V0,w,t,t0,l,Va,Dd=sp.symbols('e m d V0 omega t t0 l V_a D',positive=True)
v0=sp.sqrt(2*e*Va/m); tau=l/v0
A=e*V0/(m*d)
s=sp.symbols('s')
vx=sp.integrate(A*sp.sin(w*s),(s,t0,t))
vexit=sp.simplify(vx.subs(t,t0+tau))
xexit=sp.simplify(sp.integrate(vx,(t,t0,t0+tau)))
print('v_exit - formula:',sp.simplify(sp.expand_trig(vexit-2*A/w*sp.sin(w*tau/2)*sp.sin(w*(t0+tau/2)))))
X=xexit+vexit*Dd/v0
# D>>l approx length
Lapprox=2*2*A/w*sp.Abs(sp.sin(w*tau/2))*Dd/v0
L0=V0*l*Dd/(d*Va)
print('approx length / (L0*|sin(u)|/u):',sp.simplify((Lapprox/(L0*sp.Abs(sp.sin(w*tau/2))/(w*tau/2)))))
# at w=2pi/tau
w1=2*sp.pi/tau
Xw1=sp.simplify(X.subs(w,w1)); print('X at w1:',sp.simplify(Xw1))
amp=sp.simplify(sp.Abs(sp.diff(Xw1,t0).subs(t0,0))/w1)  # amplitude of cos/sin form
print('exact line length at w1 =',sp.simplify(2*e*V0*tau/(m*d*w1)), '=', sp.simplify((2*e*V0*tau/(m*d*w1))))
print('   in Va form:',sp.simplify((2*e*V0*tau/(m*d*w1)).subs(e,Va*0+e)), sp.simplify(sp.simplify(2*e*V0*tau/(m*d*w1))/(V0*l**2/(2*sp.pi*d*Va))))
# DC exact length
XDC=sp.simplify(A*tau**2/2+A*tau*Dd/v0); print('DC length/ (V0 l (l+2D)/(2 d Va)) =',sp.simplify(2*XDC/(V0*l*(l+2*Dd)/(2*d*Va))))
# e/m from f1
f1=sp.symbols('f1',positive=True); em=sp.solve(sp.Eq(v0, f1*l),e)[0]/m; print('e/m =',sp.simplify(em))
print('numeric e/m:', (9.4e8*0.020)**2/(2*1.0e3))
# numeric simulation of electrons with given parameters
ev,mv=1.602e-19,9.109e-31; Vav=1.0e3; lv=0.020; dv=0.010; V0v=20.0; Dv=0.30
v0v=np.sqrt(2*ev*Vav/mv); tauv=lv/v0v; print('v0',v0v,'tau',tauv,'f1',1/tauv)
def screen(wv,t0v):
    sol=solve_ivp(lambda tt,y:[y[1],ev*V0v/(mv*dv)*np.sin(wv*tt)],[t0v,t0v+tauv],[0,0],rtol=1e-12,atol=1e-20)
    xe,ve=sol.y[0,-1],sol.y[1,-1]; return xe+ve*Dv/v0v, np.max(np.abs(sol.y[0]))
for wv in [2*np.pi*1e6, np.pi/tauv, 2*np.pi/tauv]:
    T=2*np.pi/wv; t0s=np.linspace(0,T,721)
    res=[screen(wv,tt) for tt in t0s]; Xs=[r[0] for r in res]
    Lnum=max(Xs)-min(Xs)
    sub={e:ev,m:mv,d:dv,V0:V0v,w:wv,l:lv,Va:Vav,Dd:Dv}
    print('w*tau=%.3f'%(wv*tauv),' Lnum(mm)=%.5f'%(Lnum*1e3),' approx(mm)=%.5f'%(float(Lapprox.subs(sub))*1e3),' max|x| inside(mm)=%.4f'%(max(r[1] for r in res)*1e3))
print('exact at w1 formula (mm):', V0v*lv**2/(2*np.pi*dv*Vav)*1e3, ' DC exact (mm):',V0v*lv*(lv+2*Dv)/(2*dv*Vav)*1e3)
