# S0113-1 回転ドラム内壁の小物体：滑り出しと離れの順序，最も滑りやすい角
import sympy as sp, numpy as np
th,mu,k,g,m=sp.symbols('theta mu k g m',positive=True)  # k = R w^2 / g
N=m*g*(k+sp.cos(th)); f=m*g*sp.sin(th)
# (2) N=0 の角で必要な摩擦 f>0 -> mu N=0 < f : 先に滑る
thd=sp.acos(-k)
print('f at N=0:',sp.simplify(f.subs(th,thd)))  # m g sqrt(1-k^2) >0 for k<1
# (3) 条件 sin - mu cos <= mu k  for all th : max = sqrt(1+mu^2)
print('max of sin-mu cos:',sp.sqrt(1+mu**2),'at th = pi/2+atan(mu)')
x=np.linspace(0,2*np.pi,200001)
for muv in [0.3,0.5,1.0]:
    h=np.sin(x)-muv*np.cos(x); i=np.argmax(h)
    assert abs(x[i]-(np.pi/2+np.arctan(muv)))<1e-4 and abs(h[i]-np.sqrt(1+muv**2))<1e-8
# (4) mu=0.5, k=1.5 -> first angle where sin-0.5cos = 0.75
muv,kv=0.5,1.5
h=np.sin(x)-muv*np.cos(x)-muv*kv
i=np.argmax(h>0); print('slip angle deg (numeric)',np.degrees(x[i]))
th_exact=np.arcsin(muv*kv/np.sqrt(1+muv**2))+np.arctan(muv)
print('exact',np.degrees(th_exact)); assert abs(x[i]-th_exact)<1e-4
print('critical k for mu=0.5:',np.sqrt(1+muv**2)/muv)
# N>0 at slip angle
print('N/mg at slip',kv+np.cos(th_exact))
# 一般: k<1 のとき滑り角 < 離れる角 を数値で多数確認
for muv in [0.1,0.5,1,3]:
    for kv in [0.2,0.5,0.9]:
        hh=np.sin(x)-muv*(np.cos(x)+kv); i=np.argmax(hh>0)
        assert x[i]<np.arccos(-kv)
print('ok')
