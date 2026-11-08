# O20062-1 : 一様磁場中の荷電粒子にソレノイドの磁束変化で力積を与える
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
q,m,B0,b,r0,dB,R=sp.symbols('q m B0 b r0 DeltaB R',positive=True)
# (1) 誘導電場の力積: 起電力(反時計回り正) = -pi b^2 * dB1/dt ; E_theta = V/(2 pi r)
delta=sp.symbols('delta')  # 外部磁場(ソレノイド内部の追加分)の変化量(符号付き)
J_theta=lambda r,dl: -q*b**2*dl/(2*r)   # theta方向(反時計回り)の力積
v0=sp.simplify(J_theta(r0,-dB)/m); print('(1) v0 =',v0)
rho=sp.simplify(m*v0/(q*B0)); print('(2) rho=',rho,' R=r0+2rho=',sp.expand(r0+2*rho))
Rq=r0+2*rho
# (3) Qでは速度 -y (= -theta方向), 大きさ v0. 止めるには theta方向力積 +m v0
d3=sp.solve(sp.Eq(J_theta(Rq,delta),m*v0),delta)[0]
print('(3) delta =',sp.factor(d3),' ; -dB*R/r0 =',sp.simplify(d3+dB*Rq/r0))
# (4) 増加させた場合 rho 同じ, 中心 (r0-rho,0)
# 数値検証
def sim(qv,mv,B0v,bv,r0v,deltas,tswitch,T,dt=1e-4):
    # deltas: list of (t_start, dB_signed) ; ramp over dt
    def Eth(r,t):
        s=0.0
        for (ts,dl) in deltas:
            if ts<=t<ts+dt: s+=dl/dt
        return -bv**2*s/(2*r)
    def f(t,y):
        x,yy,vx,vy=y; r=np.hypot(x,yy); e=Eth(r,t)
        Ex=-e*yy/r; Ey=e*x/r
        ax=qv/mv*(Ex+vy*B0v); ay=qv/mv*(Ey-vx*B0v)
        return [vx,vy,ax,ay]
    ts=sorted(set([0]+[d[0] for d in deltas]+[d[0]+dt for d in deltas]+[T]))
    y=[r0v,0,0,0]; traj=[]
    for a,c in zip(ts[:-1],ts[1:]):
        s=solve_ivp(f,(a,c),y,rtol=1e-10,atol=1e-12,max_step=(c-a)/50)
        y=s.y[:,-1]; traj.append(s)
    return y,traj
qv,mv,B0v,bv,r0v,dBv=1.0,1.0,1.0,0.5,1.0,2.0
v0n=qv*bv**2*dBv/(2*mv*r0v); rhon=mv*v0n/(qv*B0v); Rn=r0v+2*rhon
Th=np.pi*mv/(qv*B0v)
y,_=sim(qv,mv,B0v,bv,r0v,[(0,-dBv)],None,1e-3)
print('num v after kick',y[2:],' theory (0,',v0n,')')
y,_=sim(qv,mv,B0v,bv,r0v,[(0,-dBv)],None,Th+1e-4/2)
print('num pos at half period',y[:2],' theory R=',Rn)
d3n=float(d3.subs({q:qv,m:mv,B0:B0v,b:bv,r0:r0v,dB:dBv}))
y,_=sim(qv,mv,B0v,bv,r0v,[(0,-dBv),(Th,d3n)],None,Th+1e-4+0.5)
print('after second kick, velocity',y[2:],' pos',y[:2],' delta3=',d3n,' -dB R/r0=',-dBv*Rn/r0v)
# (4) 判定: 増加 dB で最近接距離
for dBv4 in [0.5,1.5,2.5,3.5,6.0,7.0]:
    v=qv*bv**2*dBv4/(2*mv*r0v); rh=mv*v/(qv*B0v)
    # simulate full circle, find min r
    def f(t,y):
        x,yy,vx,vy=y; return [vx,vy,qv/mv*vy*B0v,-qv/mv*vx*B0v]
    s=solve_ivp(f,(0,2*Th),[r0v,0,0,-v],rtol=1e-10,atol=1e-12,max_step=0.001)
    rmin=np.min(np.hypot(s.y[0],s.y[1]))
    safe= dBv4<r0v*(r0v-bv)*B0v/bv**2 or dBv4>r0v*(r0v+bv)*B0v/bv**2
    print('dB=',dBv4,' rmin=',round(rmin,4),' >b?',rmin>bv,' theory safe:',safe)
