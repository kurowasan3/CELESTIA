# PHY-PHO-016 検証：コンプトン散乱線のブラッグ分光器での角度分離
import sympy as sp, numpy as np
n,d,lam,th,dl=sp.symbols('n d lambda theta Delta_lambda',positive=True)
# (1) n lam = 2 d sin th を微分: n dl = 2 d cos th dth -> dth = dl/lam * tan th
dth=n*dl/(2*d*sp.cos(th))
print('(1)',sp.simplify(dth.subs(n,2*d*sp.sin(th)/lam)-dl/lam*sp.tan(th)))
lam0=71.1e-12; lc=2.43e-12; D=0.201e-9
def lamp(phi): return lam0+lc*(1-np.cos(phi))
# (2) phi=90, n=1
t0=np.arcsin(lam0/(2*D)); 
approx=lc/lam0*np.tan(t0)
t1=np.arcsin(lamp(np.pi/2)/(2*D))
print('(2) theta0 deg',np.degrees(t0),' approx dtheta deg',np.degrees(approx),' exact',np.degrees(t1-t0))
# (3) theta<=75deg
for k in range(1,7):
    s0=k*lam0/(2*D); s1=k*lamp(np.pi/2)/(2*D)
    if s1<1:
        a0=np.degrees(np.arcsin(s0)); a1=np.degrees(np.arcsin(s1)); print('  n',k,a0,a1,a1-a0, 'ok' if a1<=75 else 'over')
    else: print('  n',k,'impossible',s0,s1)
t05=np.arcsin(5*lam0/(2*D)); print('(3) approx formula n=5',np.degrees(lc/lam0*np.tan(t05)))
# (4) shifted line at 64.0 deg in 5th order
lp=2*D*np.sin(np.radians(64.0))/5; dL=lp-lam0
phi=np.degrees(np.arccos(1-dL/lc)); print('(4) dlam pm',dL*1e12,' phi',phi)
# (5) visible
r=lc/500e-9; print('(5) rel',r,' dtheta at 75 deg (rad)',r*np.tan(np.radians(75)),' deg',np.degrees(r*np.tan(np.radians(75))))
