# T20253-1 気体ばねに繰り返し衝突する台車（最大圧縮時に熱浴で冷却）
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
K,U0,T0,V0,R,n=sp.symbols('K U_0 T_0 V_0 R n',positive=True)
U0e=sp.Rational(3,2)*R*T0
T1=T0*(1+K/U0e); V1=V0*(T0/T1)**sp.Rational(3,2)
T4=T0*(V1/V0)**sp.Rational(2,3)
Kp=sp.Rational(3,2)*R*(T0-T4)
print("K' =",sp.simplify(Kp), ' == U0K/(U0+K)?', sp.simplify(Kp-U0e*K/(U0e+K)))
print('1/K\'-1/K =',sp.simplify(1/Kp-1/K))
# 熱の収支
Qout=sp.Rational(3,2)*R*(T1-T0); Qin=sp.Rational(3,2)*R*(T0-T4)
print('Qout=K?',sp.simplify(Qout-K),' Qin=K\'?',sp.simplify(Qin-Kp))
# 数値：台車とピストンの運動を直接積分（R=1,T0=1,V0=1,断面積1 => 初期位置x=1）
m=1.0; Rn=1.0; T0n=1.0; V0n=1.0
def push(K,Tstart,Vstart,v0,compress):
    # 気体：p V^{5/3}=const（断熱），体積 V=x（断面積1）
    C=Rn*Tstart/Vstart*Vstart**(5/3)
    def f(t,y): x,v=y; p=C/x**(5/3); return [v,(p)/m]   # 台車は+x向き（体積増）に押される
    if compress:
        ev=lambda t,y: y[1]; ev.terminal=True; ev.direction=1
        s=solve_ivp(f,[0,100],[Vstart,-v0],events=ev,rtol=1e-11,atol=1e-13)
    else:
        ev=lambda t,y: y[0]-V0n; ev.terminal=True; ev.direction=1
        s=solve_ivp(f,[0,100],[Vstart,0.0],events=ev,rtol=1e-11,atol=1e-13)
    return s.y[:,-1]
Kn=1.5*Rn*T0n  # K0=U0
U0n=1.5*Rn*T0n
Ks=[Kn]
for i in range(5):
    v0=np.sqrt(2*Ks[-1]/m)
    x1,_=push(Ks[-1],T0n,V0n,v0,True)          # 最大圧縮体積
    T1n=T0n*(V0n/x1)**(2/3)
    # 保持して T0 に冷却（体積 x1）、解放して膨張
    x2,v2=push(None,T0n,x1,0,False)
    Ks.append(0.5*m*v2**2)
    print(f'collision {i+1}: Vmax-compr={x1:.6f} T1={T1n:.6f} (theory {T0n*(1+Ks[-2]/U0n):.6f})  K={Ks[-1]:.6f}  theory={U0n*Ks[-2]/(U0n+Ks[-2]):.6f}  K0/(n+1)={Kn/(i+2):.6f}')
