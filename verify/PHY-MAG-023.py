# Q186-1 直線電流のまわりの電子：v_z = v_z0 - K ln(r/r0)，速さ一定
import numpy as np
from scipy.integrate import solve_ivp
e=1.602e-19; me=9.109e-31; mu0=4e-7*np.pi
I=5.0; r0=2.0e-3; v0=1.0e5
K=e*mu0*I/(2*np.pi*me); print('K=',K)
def f(t,s):
    x,y,z,vx,vy,vz=s; r2=x*x+y*y
    B=mu0*I/(2*np.pi*r2)*np.array([-y,x,0.0])   # 電流 +z
    F=-e*np.cross([vx,vy,vz],B)
    return [vx,vy,vz,*(F/me)]
T=2e-7
# (A) 導線に向かって打ち出す
sA=solve_ivp(f,[0,T],[r0,0,0,-v0,0,0],rtol=1e-11,atol=1e-16,max_step=1e-10)
rA=np.hypot(sA.y[0],sA.y[1]); print('A: rmin num',rA.min(),' theory',r0*np.exp(-v0/K))
# v_z と r の関係
i=len(rA)//7; print(' vz check',sA.y[5][i],' theory',-K*np.log(rA[i]/r0))
# (B+) 電流と同じ向き
sB=solve_ivp(f,[0,T],[r0,0,0,0,0,v0],rtol=1e-11,atol=1e-16,max_step=1e-10)
rB=np.hypot(sB.y[0],sB.y[1]); print('B+: r range',rB.min(),rB.max(),' theory rmax',r0*np.exp(2*v0/K))
# (B-) 逆向き
sC=solve_ivp(f,[0,T],[r0,0,0,0,0,-v0],rtol=1e-11,atol=1e-16,max_step=1e-10)
rC=np.hypot(sC.y[0],sC.y[1]); print('B-: r range',rC.min(),rC.max(),' theory rmin',r0*np.exp(-2*v0/K))
sp=np.sqrt(sB.y[3]**2+sB.y[4]**2+sB.y[5]**2); print('speed const',sp.min(),sp.max())
print('rmax mm',r0*np.exp(2*v0/K)*1e3,' rmin mm',r0*np.exp(-2*v0/K)*1e3, 'A rmin mm', r0*np.exp(-v0/K)*1e3)
print('OK')
