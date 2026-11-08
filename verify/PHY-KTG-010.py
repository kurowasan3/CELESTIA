# P20 類題：温度の異なる2つの壁の間を往復する分子（1次元，希薄）の力
import sympy as sp, numpy as np
m,L,u1,u2,u0,N=sp.symbols('m L u1 u2 u0 N',positive=True)
T=L/u1+L/u2
F=sp.simplify(m*(u1+u2)/T); print('force per molecule (both walls)=',F)
KE=sp.simplify((sp.Rational(1,2)*m*u1**2*(L/u1)+sp.Rational(1,2)*m*u2**2*(L/u2))/T); print('time-avg KE=',KE)
# 加熱直後の右壁
Fr=sp.simplify(N*u0/(2*L)*m*(u0+u2)); print('right wall just after heating:',Fr,' vs before',N*m*u0**2/L)
# 数値シミュレーション：分子を一様に配置・向きランダム，時間平均の力
rng=np.random.default_rng(1)
def simulate(n=2000,Lv=1.0,v1=1.0,v2=3.0,v0=2.0,Ttot=200.0,t_heat=0.0):
    x=rng.uniform(0,Lv,n); v=rng.choice([-1,1],n)*v0
    impL=0.0; impR=0.0; t=0.0
    # event-driven per molecule independent: compute impulses analytically per molecule
    impR_early=0.0
    for i in range(n):
        xi,vi,ti=x[i],v[i],0.0
        while True:
            if vi>0: dt=(Lv-xi)/vi; wall='R'
            else: dt=xi/(-vi); wall='L'
            if ti+dt>Ttot: break
            ti+=dt
            if wall=='R':
                imp=m_*(abs(vi)+v2); impR+=imp; xi=Lv; vi=-v2
                if ti<Lv/v0*0.5: impR_early+=imp
            else:
                imp=m_*(abs(vi)+v1); impL+=imp; xi=0.0; vi=v1
    return impL/Ttot, impR/Ttot, impR_early/(Lv/v0*0.5)
m_=1.0
fL,fR,fRe=simulate()
print('sim force L',fL,' R',fR,' theory',2000*1.0*1.0*3.0/1.0)
print('sim right wall early',fRe,' theory',2000*2.0/(2*1.0)*1.0*(2.0+3.0))
assert abs(fL/6000-1)<0.02 and abs(fR/6000-1)<0.02 and abs(fRe/10000-1)<0.05
