# M21-1: projectile inside a truck accelerating with alpha (forward = +x)
import sympy as sp, numpy as np
from scipy.optimize import minimize_scalar
from scipy.integrate import solve_ivp
g,al,u,b=sp.symbols('g alpha u beta',positive=True)
t=2*u*sp.sin(b)/g
xr=u*sp.cos(b)*t-al*t**2/2
print('(2) vertical launch:',sp.simplify(xr.subs(b,sp.pi/2)))
# (3) launch along apparent up: tan(beta')=alpha/g tilt forward from vertical -> xr should be 0
bb=sp.pi/2-sp.atan(al/g); print('(3) xr along apparent vertical',sp.simplify(xr.subs(b,bb)))
gp=sp.sqrt(g**2+al**2)
for gv,av,uv in [(9.8,3.0,5.0),(9.8,7.0,4.0)]:
    f=sp.lambdify(b,xr.subs({g:gv,al:av,u:uv}))
    rf=minimize_scalar(lambda x:-f(x),bounds=(0,np.pi/2),method='bounded',options={'xatol':1e-12})
    rb=minimize_scalar(lambda x:f(x),bounds=(np.pi/2,np.pi),method='bounded',options={'xatol':1e-12})
    gpn=np.hypot(gv,av); phi=np.arctan(av/gv)
    print(f'alpha={av}: fwd max {-rf.fun:.6f} th {uv**2/(gpn+av):.6f}; angle {np.degrees(rf.x):.5f} th {45-np.degrees(phi)/2:.5f}')
    print(f'          bwd max {-rb.fun:.6f} th {uv**2/(gpn-av):.6f}; angle from back floor {180-np.degrees(rb.x):.5f} th {45+np.degrees(phi)/2:.5f}')
    # ground-frame integration check for the forward optimum (truck already moving with speed 3 m/s)
    V0=3.0; be=rf.x
    def rhs(tm,y): return [y[2],y[3],0,-gv]
    ev=lambda tm,y: y[1]; ev.terminal=True; ev.direction=-1
    r=solve_ivp(rhs,[0,10],[0,0,V0+uv*np.cos(be),uv*np.sin(be)],events=ev,rtol=1e-12,atol=1e-12)
    T=r.t[-1]; truck=V0*T+0.5*av*T**2
    print('   ground-frame landing rel to launch point on truck:',r.y[0,-1]-truck)
# --- (5) L_f L_b = u^4/g^2（alpha によらない）
Lf=u**2/(gp+al); Lb=u**2/(gp-al)
print('Lf*Lb =',sp.simplify(Lf*Lb)); assert sp.simplify(Lf*Lb-u**4/g**2)==0
# 一般化：床に垂直 gperp, 平行 gpar の一様な加速度場での最大到達距離の積 = u^4/gperp^2（数値最適化で確認）
for gperp,gpar,uv in [(9.8,3.0,5.0),(4.0,7.0,2.0),(1.0,0.2,1.0)]:
    psi=np.linspace(1e-6,np.pi-1e-6,400001); tt=2*uv*np.sin(psi)/gperp
    x=uv*np.cos(psi)*tt-gpar*tt**2/2
    prod=x.max()*(-x.min()); print(' gperp',gperp,'gpar',gpar,' product',prod,' u^4/gperp^2',uv**4/gperp**2)
    assert abs(prod/(uv**4/gperp**2)-1)<1e-6
print('(5) OK')
