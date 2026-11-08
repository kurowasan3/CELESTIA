# Q144-1: 片面が弾性・片面が完全非弾性の板が粒子の流れから受ける力と終端速度
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
u,v0,n,m,S,F0=sp.symbols('u v0 n m S F0',real=True)
K=n*m*S
F_mid=K*((v0-u)**2-2*(v0+u)**2)          # -v0<u<v0
F_hi=-2*K*((u+v0)**2+(u-v0)**2)           # u>v0
F_lo= K*((v0-u)**2+(u+v0)**2)             # u<-v0
print('連続性 u=v0:',sp.simplify(F_mid.subs(u,v0)-F_hi.subs(u,v0)),' u=-v0:',sp.simplify(F_mid.subs(u,-v0)-F_lo.subs(u,-v0)))
us=sp.solve(sp.Eq(F_mid,0),u); print('(2) 終端速度候補',us)
ust=[s for s in us if abs(float(s.subs(v0,1)))<1][0]; print('  ->',sp.nsimplify(ust))
assert sp.simplify(ust+(3-2*sp.sqrt(2))*v0)==0
# (3) 右向き一定の力 F0
k=sp.Symbol('k',positive=True)  # F0 = k*K*v0^2
um=sp.solve(sp.Eq(F_mid+k*K*v0**2,0),u); print('(3) 中間領域',um)
uh=sp.solve(sp.Eq(F_hi+k*K*v0**2,0),u); print('(3) 高速領域',uh)
print('切替 F0 =',sp.simplify(-F_mid.subs(u,v0)/(K*v0**2)),'K v0^2  左向きの切替 =',sp.simplify(F_mid.subs(u,-v0)/(K*v0**2)))
# 粒子シミュレーションで F(u) を測る（板の速度 u を固定，板の系で粒子を動かし衝突の力積を数える）
rng=np.random.default_rng(1)
def measure(uu,v=1.0,dens=2000.0,Lb=3.0,T=40.0,dt=2e-3):
    N=int(dens*2*Lb)
    pos=[rng.uniform(-Lb,Lb,N),rng.uniform(-Lb,Lb,N)]
    vel=[+v,-v]; imp=0.0
    for step in range(int(T/dt)):
        for c in (0,1):
            w=vel[c]-uu                       # 板に対する相対速度
            x=pos[c]; xn=x+w*dt
            if w<0:   # 右面(x>0側)へ左向きに当たる: 弾性
                hit=(x>0)&(xn<=0); imp+= -np.sum(hit)*2*(-w)   # 板が受ける力積(m=1)
            else:     # 左面(x<0側)へ右向きに当たる: 完全非弾性（粒子は板の速度になって横に去る）
                hit=(x<0)&(xn>=0); imp+= np.sum(hit)*w
            # 当たった粒子は衝突後に横へ去ったものとして数え，同じ種類の新しい粒子として
            # そのまま通過させる（周期境界で一様な密度の流れを保つ）
            xn=np.where(xn>Lb,xn-2*Lb,xn); xn=np.where(xn<-Lb,xn+2*Lb,xn)
            pos[c]=xn
    return imp/T/dens   # n=dens, m=1, S=1 で割った F/(n m S)
for uu in (-1.5,-0.5,ust.subs(v0,1),0.3,1.5):
    uu=float(uu); Fs=measure(uu)
    Ft=float((F_mid if -1<uu<1 else (F_hi if uu>1 else F_lo)).subs({u:uu,v0:1,n:1,m:1,S:1}))
    print('u=%.3f  F/(nmS) sim=%.3f  theory=%.3f'%(uu,Fs,Ft))
    assert abs(Fs-Ft)<0.03*max(1,abs(Ft))
# 運動方程式を数値積分して，力F0=3 K v0^2 と 12 K v0^2 での終端速度を確認
Ff=sp.lambdify(u,sp.Piecewise((F_hi,u>v0),(F_lo,u<-v0),(F_mid,True)).subs({v0:1,n:1,m:1,S:1}))
for kk in (0,3,12,-5):
    so=solve_ivp(lambda tt,Y:[Ff(Y[0])+kk],[0,200],[0.0],rtol=1e-10)
    print('F0=%d: u_inf=%.5f'%(kk,so.y[0][-1]))
print('理論 F0=3:',float(-3+sp.sqrt(11)),' F0=12:',float(sp.sqrt(2)),' F0=-5:', float(sp.sqrt(5/2-1))*-1)
print('OK')
