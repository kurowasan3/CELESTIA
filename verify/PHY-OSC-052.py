# S0120-1 一定の速さで進む壁のばねに繰り返し押される物体（粗い床）：周期運動と平均仕事率
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,k,g,mu,mu2,v0=sp.symbols('m k g mu muprime v0',positive=True)
x0=mu*m*g/k; x1=mu2*m*g/k
u2=v0**2+k/m*((x0-x1)**2-x1**2)
print('u^2 =',sp.expand(u2))
A=sp.sqrt(m/k*v0**2+(x0-x1)**2); print('amplitude A=',A)
# 数値シミュレーション（床の系，左向き正）
def simulate(mv,kv,gv,muv,mu2v,v0v,tend):
    # 状態: P の位置 X (左向き正), 速度 V；壁のばね自由端の位置 S(t)=s0+v0 t
    # 接触: 圧縮 c = S - X > 0 のとき力 k c（左向き）
    s0=-0.0   # t=0 で自由端が P に接する
    t=0.0; X=0.0; V=0.0; work=0.0; out=[]
    dt=2e-5
    # 半陰的オイラーで静止摩擦を正確に扱う
    n=int(tend/dt); Ws=0.0
    for i in range(n):
        S=s0+v0v*t; c=S-X; Fs=kv*c if c>0 else 0.0
        if abs(V)<1e-12:
            if Fs<=muv*mv*gv: a=0.0; V=0.0
            else: a=(Fs-mu2v*mv*gv)/mv
        else:
            a=(Fs-np.sign(V)*mu2v*mv*gv)/mv
        Vn=V+a*dt
        if V>0 and Vn<0 and Fs<=muv*mv*gv: Vn=0.0
        X+=Vn*dt; V=Vn; Ws+=Fs*v0v*dt; t+=dt
        out.append((t,X,V,c))
    return np.array(out),Ws
mv,kv,gv,muv,mu2v,v0v=1.0,40.0,9.8,0.6,0.2,0.5
arr,Ws=simulate(mv,kv,gv,muv,mu2v,v0v,12.0)
tt,X,V,c=arr.T
u=np.sqrt(v0v**2+muv*(muv-2*mu2v)*mv*gv**2/kv)
il=[i for i in range(1,len(tt)) if c[i]<=0 and c[i-1]>0][0]
print('theory leave speed',v0v+u,' numeric leave V',V[il])
assert abs(V[il]-(v0v+u))<5e-3
print('avg power numeric',Ws/tt[-1],' theory',mu2v*mv*gv*v0v)
print('avg speed numeric',X[-1]/tt[-1])
# 周期の比較（理論）
w=np.sqrt(kv/mv); X0=muv*mv*gv/kv; X1=mu2v*mv*gv/kv; Aa=np.sqrt(v0v**2/w**2+(X0-X1)**2)
psi0=np.arcsin((X0-X1)/Aa); psi1=np.pi+np.arcsin(X1/Aa)
tp=(psi1-psi0)/w; ts=(v0v+u)/(mu2v*gv); sgap=u*ts-0.5*mu2v*gv*ts**2; tw=sgap/v0v
Tc=X0/v0v+tp+ts+tw; print('cycle period theory',Tc)
# 数値で静止開始時刻を検出
starts=[tt[i] for i in range(1,len(tt)) if V[i]>0 and V[i-1]==0]
print('numeric start times',starts[:4], 'diffs',np.diff(starts[:4]))
# 長時間平均（整数周期で）
n=int(tt[-1]//Tc); iend=np.searchsorted(tt,n*Tc)
print('avg power over',n,'cycles',np.sum(kv*np.clip(c[:iend],0,None)*v0v)*2e-5/(tt[iend]),' theory',mu2v*mv*gv*v0v)
assert abs(np.sum(kv*np.clip(c[:iend],0,None)*v0v)*2e-5/(tt[iend])-mu2v*mv*gv*v0v)<0.02*mu2v*mv*gv*v0v
# mu/2<mu'<mu の場合
mu2v=0.4; v0v=1.0
arr,Ws=simulate(mv,kv,gv,muv,mu2v,v0v,15.0)
tt,X,V,c=arr.T
u1=np.sqrt(v0v**2-muv*(2*mu2v-muv)*mv*gv**2/kv)
# 2回目以降の接触時の速さ
cont=[i for i in range(1,len(tt)) if c[i]>0 and c[i-1]<=0]
print('contact speeds',[V[i] for i in cont[:4]],' theory v0-u1',v0v-u1)
assert all(abs(V[i]-(v0v-u1))<5e-3 for i in cont[1:4])
print('avg power (mu2=0.4)',Ws/tt[-1],' theory',mu2v*mv*gv*v0v)
