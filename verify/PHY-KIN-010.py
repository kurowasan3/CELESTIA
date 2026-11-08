# J004-1 検証：エレベーター内の台から飛び出す小球の斜面上の到達距離（加速度が飛び出しの瞬間に切り替わる場合）
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
H,g,a,th,v0,t=sp.symbols('H g a theta v0 t',positive=True)
def CQ(vC2,geff):
    vC=sp.sqrt(vC2); T=2*vC*sp.tan(th)/geff; return sp.simplify(vC*T/sp.cos(th))
base=CQ(2*g*H,g); print('(1)',base)
print('(2) a all along:',sp.simplify(CQ(2*(g+a)*H,g+a)))
print('(3) switch on at C:',sp.simplify(CQ(2*g*H,g+a)/base))
print('(4) switch off at C:',sp.simplify(CQ(2*(g+a)*H,g)/base))
# 着地時の速度の向き：tan(phi)=2 tan(theta)（初速によらない）
vC=sp.symbols('v_C',positive=True); ge=sp.symbols('g_e',positive=True)
T=2*vC*sp.tan(th)/ge; print('vy/vx at Q =',sp.simplify(ge*T/vC))
# (5) a を求める逆問題: 比 = 3/4 となるのは (3)型で g/(g+a)=3/4 -> a=g/3
print('(5) a:',sp.solve(sp.Eq(g/(g+a),sp.Rational(3,4)),a))
# 数値：地上系でエレベーターの運動と小球の運動を積分し，斜面との交点を求める
gn,Hn,thn,an=9.8,0.5,np.radians(30),3.0
def run(a_slide,a_fly):
    # 滑走部はエレベーター系でエネルギー保存 vC^2 = 2(g+a_slide)H（なめらかな曲面）
    vCn=np.sqrt(2*(gn+a_slide)*Hn)
    # 地上系：時刻0で C を水平に飛び出す。エレベーター速度 U0 は任意(=1.3)とし，以後加速度 a_fly
    U0=1.3
    def rhs(t,y):
        # y = [xb, yb, vxb, vyb, ye, vye]
        return [y[2],y[3],0,-gn,y[5],a_fly]
    def hit(t,y):
        # 斜面（エレベーター内）：C を原点に下向き tan(theta)
        xr=y[0]; yr=y[1]-y[4]
        return yr + np.tan(thn)*xr if t>1e-9 else 1
    hit.terminal=True; hit.direction=-1
    s=solve_ivp(rhs,[0,5],[0,0,vCn,U0,0,U0],events=hit,rtol=1e-11,atol=1e-12)
    te=s.t_events[0][0]; y=s.y_events[0][0]; xr=y[0]
    return xr/np.cos(thn)
b=run(0,0); print('numeric base',b,' formula',4*Hn*np.tan(thn)/np.cos(thn))
print('numeric (2)',run(an,an)/b,' (3)',run(0,an)/b,' formula',gn/(gn+an),' (4)',run(an,0)/b,' formula',(gn+an)/gn)
