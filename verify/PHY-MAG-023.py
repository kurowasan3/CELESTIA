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
# (3) new: drift per radial cycle, Delta z = -2 r0 K int_0^a u (e^u - e^-u)/sqrt(v0^2-K^2u^2) du < 0
from scipy.integrate import quad
a=v0/K
# substitute u = a sin(th) to remove endpoint singularity
Iu=quad(lambda th: a*np.sin(th)*(np.exp(a*np.sin(th))-np.exp(-a*np.sin(th))) /K , 0, np.pi/2)[0]   # = int u(e^u-e^-u)/sqrt(v0^2-K^2u^2) du
dz_theory=-2*r0*K*Iu
# r-integral form: 2*int_{rmin}^{rmax} v_z/sqrt(v0^2-v_z^2) dr with r=r0*exp(a sin th)
def integrand(th):
    r=r0*np.exp(a*np.sin(th)); vz=-K*np.log(r/r0)
    return vz/(v0*np.cos(th))*r*a*np.cos(th)   # sqrt(v0^2-vz^2)=v0 cos th, dr=r a cos th dth
dz_r=2*quad(integrand,-np.pi/2,np.pi/2)[0]
print('r-integral form',dz_r,' u-form',dz_theory); assert abs(dz_r/dz_theory-1)<1e-8
# simulate longer and measure z at successive minima of r (case A)
sA2=solve_ivp(f,[0,1e-6],[r0,0,0,-v0,0,0],rtol=1e-11,atol=1e-16,max_step=2e-10,dense_output=True)
tt=np.linspace(0,1e-6,2000001); Y=sA2.sol(tt); rr=np.hypot(Y[0],Y[1])
mins=[k for k in range(1,len(rr)-1) if rr[k]<rr[k-1] and rr[k]<=rr[k+1]]
dzs=np.diff(Y[2][mins])
print('Delta z per cycle: simulated',dzs[:3],' theory',dz_theory)
assert dz_theory<0 and np.allclose(dzs,dz_theory,rtol=2e-3)
# general launch (parallel, with current): drift also negative (u measured from r_c where v_z=0)
sB2=solve_ivp(f,[0,1e-6],[r0,0,0,0,0,v0],rtol=1e-11,atol=1e-16,max_step=2e-10)
sC2=solve_ivp(f,[0,1e-6],[r0,0,0,0,0,-v0],rtol=1e-11,atol=1e-16,max_step=2e-10)
print('net z after 1us: parallel(with I)',sB2.y[2,-1],' parallel(against I)',sC2.y[2,-1],' both negative')
assert sB2.y[2,-1]<0 and sC2.y[2,-1]<0
print('OK')
