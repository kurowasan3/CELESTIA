# T20243-1 掃引音（三角波状の周波数変調）による距離と速さの測定
import numpy as np, sympy as sp
f,Vs,us=sp.symbols('f V u',positive=True)
fr=f*(Vs+us)/(Vs-us)
print('(2) exact', fr, ' series', sp.series(fr,us,0,2))
V=340.0; fc=40.0e3; Df=2.0e3; T=0.10; kap=Df/T
def fs(t):
    ph=np.mod(t,2*T)
    return np.where(ph<T, fc-Df/2+kap*ph, fc+Df/2-kap*(ph-T))
def te_of(t,L,u,t0):
    # 目標位置 x(t')=L-u(t'-t0)。時刻 t に O に戻る反射波の反射時刻 tm と発射時刻 te
    tm=(t-(L+u*t0)/V)/(1-u/V)
    return tm-(L-u*(tm-t0))/V
def beat_at(L,u,t):
    h=1e-7; te=te_of(t,L,u,t); dte=(te_of(t+h,L,u,t)-te_of(t-h,L,u,t))/(2*h)
    return abs(float(fs(t))-float(fs(te))*dte)
def approx(L,u): kt=kap*2*L/V; D=2*u*fc/V; return abs(kt-D),kt+D
for L in [1.70,1.445]:
    ex=(beat_at(L,0.51,0.05),beat_at(L,0.51,0.15)); ap=approx(L,0.51)
    print('L=',L,' exact(up,down)=',ex,' approx=',ap,' rel.err',[abs(a-b)/b for a,b in zip(ex,ap)])
# 測定値から2通りの解釈
for (bu,bd) in [(80,320),(50,290)]:
    for kt,D in [((bu+bd)/2,(bd-bu)/2),((bd-bu)/2,(bu+bd)/2)]:
        print((bu,bd),' L=',kt/kap*V/2,' u=',D*V/(2*fc))
print('consistency: (1.70-1.445)/0.5 =',(1.70-1.445)/0.5)
