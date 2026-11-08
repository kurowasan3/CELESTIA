# T20131-1 : 向かい合うばね振り子2個の弾性衝突の繰り返し（重心座標と相対座標）
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
w=1.0; d=1.0
def simulate(u,tmax,ncol=8):
    # x1,x2: 各小球の自然長位置からの変位（右向き正）. 接触: x1-x2=d
    t0=0.0; s=np.array([0,0,u,0.0]); cols=[]
    while t0<tmax and len(cols)<ncol:
        ev=lambda t,s:(s[0]-s[1])-d; ev.terminal=True; ev.direction=1
        sol=solve_ivp(lambda t,s:[s[2],s[3],-w**2*s[0],-w**2*s[1]],[t0,tmax],s,events=ev,rtol=1e-12,atol=1e-12,max_step=0.01)
        if sol.status!=1: break
        t0=sol.t_events[0][0]; s=sol.y_events[0][0].copy()
        cols.append((t0,s[0],s[1],s[2],s[3]))
        s[2],s[3]=s[3],s[2]      # 等質量の弾性衝突＝速度の交換
        s[0]-=1e-12
    return cols,s,t0
for u in [1.2,2/np.sqrt(3),1.6,3.0]:
    cols,_,_=simulate(u,200)
    ts=[c[0] for c in cols]; dts=np.diff(ts)
    s_=u/w; phi=np.arccos(d/s_)
    print('u=%.4f t1 num %.6f th %.6f | intervals'%(u,ts[0],np.arcsin(d/s_)/w),np.round(dts,6),' th',round(2*(np.pi-phi)/w,6))
u=2/np.sqrt(3)
cols,_,_=simulate(u,200,ncol=8)
for k,c in enumerate(cols[:7]):
    print('collision',k+1,'t=%.5f x1=%.5f x2=%.5f v1=%.5f v2=%.5f'%c)
# 状態の周期性: t=10π で初期状態 (0,0,u,0) に戻るか
cols,s,t0=simulate(u,10*np.pi-1e-9,ncol=100)
sol=solve_ivp(lambda t,s:[s[2],s[3],-s[0],-s[1]],[t0,10*np.pi],s,rtol=1e-12,atol=1e-12)
print('state at 10pi:',np.round(sol.y[:,-1],6),' #collisions',len(cols))
# より短い時間で戻らないこと: 2π,4π,...,8π と 5π/3 の倍数を確認
for T in [2*np.pi,4*np.pi,5*np.pi/3*3, 6*np.pi,8*np.pi]:
    cols,s,t0=simulate(u,T-1e-9,ncol=100)
    sol=solve_ivp(lambda t,s:[s[2],s[3],-s[0],-s[1]],[t0,T],s,rtol=1e-12,atol=1e-12)
    print('T=%.4f state'%T,np.round(sol.y[:,-1],4))
# 記号的: 2回目の衝突位置 x1=d cos2φ, cosφ=d/s
s_,dd=sp.symbols('s d',positive=True)
phi=sp.acos(dd/s_)
print('x1 at 2nd:',sp.simplify(dd*sp.cos(2*phi)))
