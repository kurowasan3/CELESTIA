# Q095-1 ローラー列上の板：終端速度・熱の割合・糸を切った後の減速距離
import sympy as sp, numpy as np
F,d,m,M,l,v=sp.symbols('F d m M ell v',positive=True)
# 定常：1区間 d あたりの力積 F*d/v = (m/2) v
vinf=sp.solve(sp.Eq(F*d/v,m*v/2),v)[0]; print('vinf',vinf)
# 減速：有効質量 M' dv/dx = -(m/2) v / d
Mp=M+(l/d)*m/2
x=sp.symbols('x',positive=True)
xhalf=sp.simplify(sp.log(2)*2*d*Mp/m); print('x_half',sp.expand(xhalf))
# ---- 離散シミュレーション（瞬間的な回転合わせ＝完全非弾性の結合） ----
def sim(Fv,dv,mv,Mv,lv,cut_after=None,N=200000):
    n_under=int(round(lv/dv)); Mpv=Mv+n_under*mv/2
    vv=0.0; pos=0.0; t=0.0; heat=0.0; spin=0.0; work=0.0
    hist=[]
    for i in range(N):
        Fi=Fv if (cut_after is None or i<cut_after) else 0.0
        a=Fi/Mpv
        # 次のローラーまで距離 dv 進む
        if a>0:
            tt=(-vv+np.sqrt(vv*vv+2*a*dv))/a
        else:
            tt=dv/vv
        work+=Fi*dv
        vv=vv+a*tt; t+=tt; pos+=dv
        # 新ローラー（有効質量 m/2）と結合
        vn=Mpv*vv/(Mpv+mv/2)
        heat+=0.5*Mpv*vv**2-0.5*(Mpv+mv/2)*vn**2
        spin+=0.25*mv*vn**2   # 後ろから離れるローラーが持ち去る分（定常で同数）
        vv=vn
        hist.append((pos,vv,t,work,heat,spin))
    return np.array(hist)
Fv,dv,mv,Mv,lv=2.0,0.05,0.02,5.0,2.0
h=sim(Fv,dv,mv,Mv,lv,N=40000)
vth=float(vinf.subs({F:Fv,d:dv,m:mv}))
seg=h[-2000:]
vavg=(seg[-1,0]-seg[0,0])/(seg[-1,2]-seg[0,2])
print('vinf theory',vth,' sim avg',vavg)
dW=seg[-1,3]-seg[0,3]; dH=seg[-1,4]-seg[0,4]; dS=seg[-1,5]-seg[0,5]
print('heat fraction',dH/dW,' spin fraction',dS/dW)
# 糸を切る
h2=sim(Fv,dv,mv,Mv,lv,cut_after=40000,N=80000)
pc=h2[39999,0]; vc=h2[39999,1]
idx=np.argmax(h2[40000:,1]<vc/2)+40000
print('x_half sim',h2[idx,0]-pc,' theory',float(xhalf.subs({d:dv,m:mv,M:Mv,l:lv})))
