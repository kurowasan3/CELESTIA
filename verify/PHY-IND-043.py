# Q182-1 超伝導リング：磁束保存→線形復元力→単振動
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
z,h,m,g,L,a,B0,al,I=sp.symbols('z h m g L a B0 alpha I',real=True)
Bz=B0*(1-al*z); Br=al*B0*a/2
# 磁束保存（z=h で電流0として冷却）
Isol=sp.solve(sp.Eq(sp.pi*a**2*Bz+L*I, sp.pi*a**2*B0*(1-al*h)),I)[0]
print('I(z)=',sp.simplify(Isol))
Fz=sp.simplify(-Isol*Br*2*sp.pi*a)   # dF = I dl phi x B_r r = -I B_r dl z
print('Fz=',sp.factor(Fz))
k=(sp.pi*al*a**2*B0)**2/L
assert sp.simplify(Fz-(-k*(z-h)))==0
# エネルギー法：U(z)=mgz+ (Phi0-pi a^2 Bz)^2/(2L)  -> F_mag=-d/dz(LI^2/2) at const flux
Um=L*Isol**2/2
assert sp.simplify(-sp.diff(Um,z)-Fz)==0
# h=0 の場合
zc=-m*g/k; print('最下点', sp.simplify(2*zc), ' 周期', sp.simplify(2*sp.pi*sp.sqrt(m/k)))
Imax=Isol.subs({h:0,z:2*zc}); print('Imax',sp.simplify(Imax))
# 最下点で 1/2 L I^2 = m g * 2mg/k
assert sp.simplify(L*Imax**2/2 - m*g*2*m*g/k)==0
# 浮いてつり合う条件 F(0)=mg
print('h=',sp.solve(sp.Eq(Fz.subs(z,0),m*g),h))
# 数値積分（電流の式を使わず、dPhi/dt=0 の回路方程式 L dI/dt = -pi a^2 dBz/dt を同時に解く）
vals=dict(m=5e-5,g=9.8,L=1.3e-8,a=5e-3,B0=0.5,alpha=2.0)
kk=float(k.subs({L:vals['L'],al:vals['alpha'],a:vals['a'],B0:vals['B0']}))
def f(t,y):
    zz,vv,ii=y
    F=-ii*vals['alpha']*vals['B0']*vals['a']/2*2*np.pi*vals['a']
    dI=-np.pi*vals['a']**2*(-vals['B0']*vals['alpha']*vv)/vals['L']
    return [vv,F/vals['m']-vals['g'],dI]
T=2*np.pi*np.sqrt(vals['m']/kk)
s=solve_ivp(f,[0,2*T],[0,0,0],rtol=1e-10,atol=1e-14,max_step=T/400)
print('k=',kk,'T=',T,'zmin num',s.y[0].min(),'theory',-2*vals['m']*vals['g']/kk)
assert abs(s.y[0].min()/(-2*vals['m']*vals['g']/kk)-1)<1e-4
i=np.argmin(s.y[0][:len(s.t)//2]); print('t at min',s.t[i],'T/2',T/2)
print('OK')
