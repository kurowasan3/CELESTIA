# Q077-1: 自然長で張った2本のばねの中央のおもりの横振動（3乗の復元力）と鉛直面内のたるみの振動
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
k,l,y,m,g=sp.symbols('k ell y m g',positive=True)
Fexact=2*k*(sp.sqrt(l**2+y**2)-l)*y/sp.sqrt(l**2+y**2)
print('series:',sp.series(Fexact,y,0,6))
assert sp.simplify(sp.series(Fexact,y,0,4).removeO()-k*y**3/l**2)==0
# 相似則：周期 ∝ 1/A を数値で（厳密なばね力）
def period(A,kv=50.0,lv=0.5,mv=0.1,grav=0.0,y0=0.0):
    def f(t,s):
        yy=s[0]; L=np.sqrt(lv**2+yy**2)
        F=-2*kv*(L-lv)*yy/L - mv*grav
        return [s[1],F/mv]
    ev=lambda t,s: s[1]; ev.direction=-1
    sol=solve_ivp(f,[0,200],[y0+A,0],events=ev,rtol=1e-11,atol=1e-13,max_step=0.01)
    # 速度が+から-に変わる（上端に戻る）時刻の間隔
    te=sol.t_events[0]; return te[1]-te[0] if len(te)>1 else te[0]
T1=period(0.005); T2=period(0.010)
print('T(5mm)=%.5f T(10mm)=%.5f ratio=%.5f'%(T1,T2,T1/T2))
assert abs(T1/T2-2)<0.01
# 純粋な3乗力の周期の係数：T = C * l/(A) * sqrt(m/k)
C=sp.Integral(1/sp.sqrt(1-sp.Symbol('u')**4),(sp.Symbol('u'),0,1)).evalf()
print('4*sqrt2*int =',float(4*sp.sqrt(2)*C), ' num T*A/(l sqrt(m/k))=', T1*0.005/(0.5*np.sqrt(0.1/50)))
# 鉛直面内：つり合いのたるみ y0 と微小振動の周期 2π sqrt(y0/(3g))
y0=sp.symbols('y0',positive=True)
keq=sp.diff(k*y**3/l**2,y).subs(y,y0)
w2=sp.simplify((keq/m).subs(k,m*g*l**2/y0**3)); print('omega^2 =',w2); assert sp.simplify(w2-3*g/y0)==0
kv,lv,mv,gv=200.0,0.6,0.05,9.8
# 厳密なつり合い
from scipy.optimize import brentq
yeq=brentq(lambda yy: 2*kv*(np.sqrt(lv**2+yy**2)-lv)*yy/np.sqrt(lv**2+yy**2)-mv*gv,1e-6,0.3)
yapp=(mv*gv*lv**2/kv)**(1/3)
def periodV(A):
    def f(t,s):
        yy=s[0]; L=np.sqrt(lv**2+yy**2)
        return [s[1],(-2*kv*(L-lv)*yy/L - mv*gv)/mv]
    ev=lambda t,s: s[1]; ev.direction=1
    sol=solve_ivp(f,[0,5],[-yeq-A,0],events=ev,rtol=1e-11,atol=1e-13,max_step=1e-3)
    te=sol.t_events[0]; return te[1]-te[0]
Tn=periodV(1e-5); Tth=2*np.pi*np.sqrt(yeq/(3*gv))
print('y_eq=%.5f (approx %.5f)  T_num=%.5f T_th=%.5f'%(yeq,yapp,Tn,Tth)); assert abs(Tn/Tth-1)<0.01
print('numeric example y0=0.020 m:',2*np.pi*np.sqrt(0.020/(3*9.8)))
print('OK')
