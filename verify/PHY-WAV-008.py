# P44-1 つるした重いロープを伝わる横波パルス
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp, quad
y,L,g,m,t=sp.symbols('y L g m t',positive=True)   # m = M/rho
# 次元解析 v = T^a rho^b : [N]=kg m s^-2, [rho]=kg/m
a_,b_=sp.symbols('a b'); sol=sp.solve([a_+b_, a_-b_-1+0, ],[a_,b_])  # kg: a+b=0, m: a-b=1, s: -2a=-1
print('exponents',sp.solve([a_+b_, a_-b_-1, -2*a_+1],[a_,b_]))
v=sp.sqrt(g*y)
T0=sp.integrate(1/v,(y,0,L)); print('t bottom->top',sp.simplify(T0))
# 下端から: dy/dt = sqrt(g y) -> y = g t^2/4
print('ode check', sp.simplify(sp.diff(g*t**2/4,t)-sp.sqrt(g*g*t**2/4)))
# 上端からのパルス sqrt(y)=sqrt(L)-sqrt(g)t/2, 出会い
tm=sp.solve(sp.Eq(sp.sqrt(g)*t/2, sp.sqrt(L)-sp.sqrt(g)*t/2),t)[0]
print('meet t',tm,' y',sp.simplify((g*tm**2/4)))
# おもり付き
vM=sp.sqrt(g*(y+m)); TM=sp.simplify(sp.integrate(1/vM,(y,0,L))); print('with M',TM)
mm=sp.solve(sp.Eq(TM, T0/2),m); print('half time m=',mm)
print('limit', sp.limit(TM*sp.sqrt(m),m,sp.oo))
# 数値：L=10,g=9.8
Lv,gv=10.0,9.8
print('numeric T0',quad(lambda yy:1/np.sqrt(gv*yy),0,Lv)[0], 2*np.sqrt(Lv/gv))
# 離散化した鎖（N個の質点）で横波パルスの到達を数値シミュレーション
N=400; dy=Lv/N; rho=1.0
ypos=(np.arange(N)+0.5)*dy          # 質点の高さ（下端から）
Tseg=rho*gv*np.arange(1,N)*dy        # 質点 i と i+1 の間の張力（その下の質量×g）
def rhs(tt,s):
    x=s[:N]; vx=s[N:]
    f=np.zeros(N)
    dx=np.diff(np.append(x,0.0))     # 上端（i=N）は固定 x=0
    Tall=np.append(Tseg, rho*gv*N*dy)
    force=Tall*dx/dy
    f+=force; f[1:]-=force[:-1]
    return np.concatenate([vx,f/(rho*dy)])
# 上端付近に狭いパルス（下向き）を置き，下端到達時刻を見る -> 2 sqrt(L/g) と比較
x0=np.exp(-((ypos-9.5)/0.12)**2)*1e-3
v0=np.zeros(N)
# 下向きに進む成分だけにする: 速度 = +v dx/dy (下向き進行)
vloc=np.sqrt(gv*ypos); v0=vloc*np.gradient(x0,dy)
sol=solve_ivp(rhs,(0,2.2),np.concatenate([x0,v0]),max_step=2e-3,t_eval=np.linspace(0,2.2,111))
peaks=[(tt,ypos[np.argmax(np.abs(sol.y[:N,i]))]) for i,tt in enumerate(sol.t)]
# 解析: 上端付近 y0=9.5 から下向き: sqrt(y)=sqrt(9.5)-sqrt(g) t/2
for tt,yp in peaks[::10]:
    an=(np.sqrt(9.5)-np.sqrt(gv)*tt/2)**2 if np.sqrt(9.5)-np.sqrt(gv)*tt/2>0 else 0
    print('t=%.2f sim peak y=%.2f analytic %.2f'%(tt,yp,an))
