# Q073-1: 2台のはかり（水槽とコップ）の読みの時間変化
import sympy as sp, numpy as np
g,h,v0,mu,T,MA,MB,t=sp.symbols('g h v0 mu T M_A M_B t',positive=True)
v1=sp.sqrt(v0**2+2*g*h); t1=(v1-v0)/g
# 落下時間：h = v0 t + g t^2/2
assert sp.simplify(v0*t1+g*t1**2/2-h)==0
A1=(MA-mu*t)*g-mu*v0; A2=(MA-mu*T)*g
B0=MB*g; B1=(MB+mu*(t-t1))*g+mu*v1; B2=(MB+mu*T)*g
W=(MA+MB)*g
dev1=sp.simplify(A1+B0-W); dev2=sp.simplify(A1+B1-W); dev3=sp.simplify(A2+B1-W)
print('0<t<t1:',dev1,' t1<t<T:',dev2,' T<t<T+t1:',sp.simplify(dev3))
assert sp.simplify(dev1+mu*(v0+g*t))==0 and dev2==0 and sp.simplify(dev3-mu*(v0+g*(t-T)))==0
I=sp.integrate(dev1,(t,0,t1))+sp.integrate(dev3,(t,T,T+t1))
print('time integral of deviation:',sp.simplify(I)); assert sp.simplify(I)==0
# 数値：空中の水の運動量 P(t) を水の小片で追い，合計の読み = W - dP/dt を確認
gv,hv,v0v,muv,Tv=9.8,0.80,1.0,0.30,3.0
dt=1e-4; ts=np.arange(0,Tv+2,dt)
P=[];
t1v=float(t1.subs({g:gv,h:hv,v0:v0v}))
for tt in ts:
    # 時刻 tt に空中にある小片：放出時刻 s∈[max(0,tt-t1), min(tt,T)]
    a=max(0,tt-t1v); b=min(tt,Tv)
    if b<=a: P.append(0.0); continue
    # 運動量 = ∫ mu*(v0+g(tt-s)) ds
    P.append(muv*(v0v*(b-a)+gv*((tt-a)**2-(tt-b)**2)/2))
P=np.array(P); dP=np.gradient(P,dt)
for tt in [0.1,0.2,1.0,3.1,3.3]:
    i=int(tt/dt); th=float((dev1 if tt<t1v else (dev2 if tt<Tv else dev3)).subs({g:gv,h:hv,v0:v0v,mu:muv,T:Tv,t:tt}))
    print(tt, -dP[i], th); assert abs(-dP[i]-th)<1e-3
print('t1=',t1v,'OK')
