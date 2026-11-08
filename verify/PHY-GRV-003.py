# P17-1: ホーマン遷移の出発時の惑星配置と帰りの待ち時間
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
C=sp.symbols('C',positive=True)
tH=sp.Rational(1,2)*((1+C)/2)**sp.Rational(3,2)        # 年
phi=sp.pi*(1-((1+C)/(2*C))**sp.Rational(3,2))            # 出発時に火星が先行する角
wE=2*sp.pi; wM=2*sp.pi/C**sp.Rational(3,2)
vA=sp.sqrt(2*C/(1+C))                                    # v_E 単位
print('vA/vE',vA, ' vB/vE', vA/C)
Cn=sp.Rational(152,100)
print('tH',float(tH.subs(C,Cn)),'phi deg',float(sp.deg(phi.subs(C,Cn))))
x=(2*sp.pi-2*wE*tH)
xm=float(x.subs(C,Cn))%(2*np.pi)
wait=xm/float((wE-wM).subs(C,Cn)); print('wait',wait,'total',2*float(tH.subs(C,Cn))+wait)
# 待ち時間 0 の条件
print('C0 =',sp.nsimplify(2**sp.Rational(5,3)-1), float(2**sp.Rational(5,3)-1), ' tH there', float(tH.subs(C,2**sp.Rational(5,3)-1)))
# 数値軌道計算で照合（GM=4pi^2 AU^3/yr^2）
GM=4*np.pi**2; c=1.52
def f(t,y):
    r=np.hypot(y[0],y[1]); return [y[2],y[3],-GM*y[0]/r**3,-GM*y[1]/r**3]
vE=2*np.pi; vAn=vE*np.sqrt(2*c/(1+c))
ev=lambda t,y:y[1]; ev.direction=-1; ev.terminal=True
s=solve_ivp(f,[1e-9,3],[1,0,0,vAn],events=ev,rtol=1e-12,atol=1e-12,first_step=1e-6)
tB=s.t[-1]; print('transfer time num',tB,'r at B',-s.y[0,-1])
ph=float(phi.subs(C,Cn)); thM=ph+2*np.pi/c**1.5*tB; print('Mars angle at arrival (deg, should be 180)',np.degrees(thM))
# 帰り: 到着 + wait に出発し，tH 後に地球が到着点にいるか
tdep=tB+wait; thMdep=ph+2*np.pi/c**1.5*tdep; thEarr=2*np.pi*(tdep+tB)
print('return mismatch (deg)',(np.degrees(thEarr-(thMdep+np.pi))+180)%360-180)
