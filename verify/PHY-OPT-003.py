# M75-1 検証：歯車による往復光の時間平均強度（三角波），整数の不定性の解消
import numpy as np, sympy as sp
c=3.0e8; N=200
def I_ratio(nu,l_eff,samples=400000):
    # 細い光線：外向き通過時刻 t に歯車が開いていて，t+tau にも開いていれば検出
    tau=2*l_eff/c
    Tg=1/(N*nu)                           # 歯1周期
    t=np.random.default_rng(0).uniform(0,Tg*50,samples)
    open_=lambda tt: ((tt/Tg)%1.0)<0.5     # 歯と隙間は同幅
    both=open_(t)&open_(t+tau)
    slow=0.5
    return both.mean()/slow
# (1) 理論式 |1-2s|, s=frac(2lN nu/c)
l=750.0
for nu in [100,333.3,500,1000,1600,1900]:
    K=2*l*N/c; s=(K*nu)%1
    print(nu,'sim',round(I_ratio(nu,l),3),'theory',round(abs(1-2*s),3))
# (2) データからの逆算
nuP,nuQ,IP,IQ=1600.,1900.,0.20,0.80
sP=(1+IP)/2; sQ=(1+IQ)/2          # 増加中 => s in (1/2,1)
K=(sQ-sP)/(nuQ-nuP); m=K*nuP-sP
print('K=',K,'m=',m,'l=',K*c/(2*N))
print('first darkness nu1=',1/(2*K),' darkness count below nuP:',int(np.floor(K*nuP+0.5)))
# (3) 管 L=300 m を液体で満たす，nuQ 固定
L=300.; C=nuQ*2*N*L/c
print('C=',C)
# 最大を1回通過後，減少して 0.70 => K nu = 2+0.15
n=1+(2.15-1.9)/C; print('n=',n)
# 数値：充填割合を0->1で変化させ，強度の推移を追う
xs=np.linspace(0,1,2001)
I=[abs(1-2*(((2*N*(l+(n-1)*L*x)/c)*nuQ)%1)) for x in xs]
I=np.array(I); peaks=np.sum((I[1:-1]>I[:-2])&(I[1:-1]>=I[2:]))
print('final I=',I[-1],' number of maxima during filling=',peaks)
print('sim final:',I_ratio(nuQ,l+(n-1)*L))
