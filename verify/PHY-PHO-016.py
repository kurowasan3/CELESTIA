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
r=lc/500e-9; print('(remark: visible light) rel',r,' dtheta at 75 deg (rad)',r*np.tan(np.radians(75)),' deg',np.degrees(r*np.tan(np.radians(75))))

# ---- new (5): wavelength condition for resolving Compton shift; atom recoil ----
lc_s,phi_s,tm,dth_s,lam_s=sp.symbols('lambda_C phi theta_max delta_theta lambda',positive=True)
cond=sp.solve(sp.Eq(lc_s*(1-sp.cos(phi_s))/lam_s*sp.tan(tm),dth_s),lam_s)[0]
print('lambda_max =',cond)
# illustrate: with delta_theta = 1e-3 rad, theta_max=75deg, phi=90deg
lmax=lc*(1-0)*np.tan(np.radians(75))/1e-3; print('lambda_max (dtheta=1e-3 rad) =',lmax,'m  (X-ray ok, visible 5e-7 m not)')
print('relative change visible:',lc/500e-9)
# Compton formula with recoiling mass M: derive and compare
h_,c_,Mr,E0,E1=sp.symbols('h c M E0 E1',positive=True)
# photon E0 -> E1 at angle phi on mass Mr at rest (relativistic): 1/E1 - 1/E0 = (1-cos)/(Mr c^2)
p0=E0/c_; p1=E1/c_
# energy: E0 + Mr c^2 = E1 + sqrt((Mr c^2)^2 + (p0^2+p1^2-2p0p1cos) c^2)
eq=sp.Eq((E0+Mr*c_**2-E1)**2,(Mr*c_**2)**2+(p0**2+p1**2-2*p0*p1*sp.cos(phi_s))*c_**2)
sol=sp.solve(eq,E1)[0]
dlam=sp.simplify(h_*c_/sol-h_*c_/E0); print('Delta lambda for recoil mass M =',dlam)
me=9.109e-31; Mc=12*1.6605e-27; print('m/M for carbon:',me/Mc,'  atomic shift (pm) at 90deg:',lc*me/Mc*1e12)
