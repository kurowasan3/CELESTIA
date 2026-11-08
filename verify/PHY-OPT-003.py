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
# (3)(a) 太い光線：各部分が G の同じ位置を往復する場合も I/Is=|1-2s|
rng=np.random.default_rng(1)
def I_thick(x,w,gamma=0.5,samples=400000):
    # 時間は歯の周期 Tg、位置は歯1つ分の幅 p を単位とする。歯は +y 向きに速さ 1 (p/Tg) で動く
    t=rng.uniform(0,50,samples); y=rng.uniform(0,w,samples)   # 光線内の一様な位置
    open_=lambda yy,tt: ((yy-tt)%1.0)<gamma
    return (open_(y,t)&open_(y,t+x)).mean()
for w in [0.0,0.3,0.7,1.0,2.4]:
    Is=I_thick(0.0,w)
    vals=[round(I_thick(x,w)/Is,3) for x in [0.1,0.25,0.5,0.8,1.5]]
    print('w=',w,'I/Is at x=0.1,.25,.5,.8,1.5 :',vals,' theory',[round(abs(1-2*(x%1)),3) for x in [0.1,0.25,0.5,0.8,1.5]])
    assert all(abs(v-abs(1-2*(x%1)))<0.01 for v,x in zip(vals,[0.1,0.25,0.5,0.8,1.5]))
# 対照：帰りの光が「瞬間の透過率」だけで決まると誤って考えた場合（T(t)T(t+x) の積）は s=1/2 で 0 にならない
def wrong_model(x,w,samples=200000):
    t=rng.uniform(0,50,samples)
    def T(tt):
        ys=np.linspace(0,w,201)[None,:]
        return (((ys-tt[:,None])%1.0)<0.5).mean(axis=1)
    return (T(t)*T(t+x)).mean()/(T(t)*T(t)).mean()
print('wrong (product-of-transmissions) model at s=1/2, w=0.4:',round(wrong_model(0.5,0.4,20000),3),'(nonzero: the claim relies on this)')
# (3)(b) すき間の割合 gamma
for g in [0.3,0.5,0.6,0.75,0.9]:
    Is=I_thick(0.0,0.0,g); r=I_thick(0.5,0.0,g)/Is
    th=0.0 if g<=0.5 else 2-1/g
    xs=np.linspace(0,1,101); mins=min(I_thick(x,0.0,g,40000) for x in xs)/Is
    print('gamma',g,'I/Is(s=1/2)=',round(r,3),'theory',round(th,3),' min over s:',round(mins,3))
    assert abs(r-th)<0.01 and abs(mins-th)<0.02
print('OK')
