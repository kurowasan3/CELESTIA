# 2004[3] 類題：突然吹き始める風・止む風と，伝わっている途中の波の振動数
import numpy as np, sympy as sp
def arrivals(f0,c,L,v,wfun_breaks,tmin,tmax):
    """wfun_breaks: list of (t_start, w) piecewise-constant wind (sorted). 波の山を発射時刻ごとに追跡し，x=L への到着時刻を返す"""
    out=[]
    for k in range(int(tmin*f0),int(tmax*f0)):
        te=k/f0; x=v*te; t=te
        # 区分ごとに進める
        bks=[b for b in wfun_breaks]
        while True:
            w=[ww for (ts,ww) in bks if ts<=t+1e-15][-1]
            nxt=[ts for (ts,ww) in bks if ts>t+1e-15]
            tn=nxt[0] if nxt else np.inf
            sp_=c+w
            ta=t+(L-x)/sp_
            if ta<=tn: out.append((te,ta)); break
            x+=sp_*(tn-t); t=tn
    return np.array(out)
f0,c,L=500.0,340.0,680.0
w=30.0
# 音源静止，t=0 に風が吹き始め t=T に止む
T=5.0
arr=arrivals(f0,c,L,0.0,[(-1e9,0.0),(0.0,w),(T,0.0)],-5,12)
ta=arr[:,1]; fobs=1/np.diff(ta); tm=ta[1:]
def fr_at(t): return fobs[np.argmin(abs(tm-t))]
print('t=0.5 :',fr_at(0.5),' theory',f0*(c+w)/c)
print('t=3   :',fr_at(3.0),' theory',f0)
print('t=T+0.5:',fr_at(T+0.5),' theory',f0*c/(c+w))
print('t=T+3 :',fr_at(T+3),' theory',f0)
# 切り替わる時刻
hi=tm[(fobs>f0*1.01)]; print('high interval',hi.min(),hi.max(),' theory 0 ..',L/(c+w))
lo=tm[(fobs<f0*0.99)]; print('low interval',lo.min(),lo.max(),' theory',T,'..',T+L/c)
# 余分な波の数
n_extra=np.sum((ta>0)&(ta<L/(c+w)))-f0*L/(c+w)
print('extra waves',n_extra,' theory',f0*L*w/(c*(c+w)))
assert abs(n_extra-f0*L*w/(c*(c+w)))<1.5
# 動く音源
v=20.0
arr=arrivals(f0,c,L,v,[(-1e9,0.0),(0.0,w)],-20,10)
ta=arr[:,1]; fobs=1/np.diff(ta); tm=ta[1:]
print('before',fobs[np.argmin(abs(tm+0.5))],f0*c/(c-v))
print('mid',fobs[np.argmin(abs(tm-0.5))],f0*(c+w)/(c-v))
print('after',fobs[np.argmin(abs(tm-(L/(c+w)+0.5)))],f0*(c+w)/(c+w-v))
# 記号での確認
f0s,cs,ws,Ls=sp.symbols('f0 c w L',positive=True)
print(sp.simplify(f0s*(cs+ws)/cs*Ls/(cs+ws)-f0s*Ls/(cs+ws)-f0s*(Ls/cs-Ls/(cs+ws))))
