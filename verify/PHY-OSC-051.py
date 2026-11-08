# S0114-1 穴を通した糸の円すい振り子をゆっくり引き上げる（軸まわりの角運動量保存と断熱不変量）
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
from scipy.optimize import brentq
l,th,g,m=sp.symbols('l theta g m',positive=True)
w=sp.sqrt(g/(l*sp.cos(th))); r=l*sp.sin(th)
L=m*r**2*w
print('L^2/(m^2 g)=',sp.simplify(L**2/(m**2*g)))
# (3) theta0=60deg, l->l0/2 : sin^4/cos = 8 * sin^4(60)/cos(60) = 9
F=lambda t: np.sin(t)**4/np.cos(t)-9
t1=brentq(F,0.5,np.pi/2-1e-9); print('theta1 deg',np.degrees(t1),'cos',np.cos(t1))
E=lambda ll,t: ll*(np.sin(t)**2/(2*np.cos(t))-np.cos(t))  # in units of m g
W=E(0.5,t1)-E(1.0,np.pi/3); print('W/(m g l0)=',W)
# 独立検証1: W = ∫ T (-dl), T = mg/cos(theta(l))
def th_of(ll):
    C=np.sin(np.pi/3)**4/np.cos(np.pi/3)  # at l=1
    return brentq(lambda t: ll**3*np.sin(t)**4/np.cos(t)-C,1e-6,np.pi/2-1e-12)
ls=np.linspace(0.5,1.0,20001); Ts=np.array([1/np.cos(th_of(x)) for x in ls])
Wint=np.trapezoid(Ts,ls); print('W by integral of tension',Wint); assert abs(Wint-W)<1e-6
# 独立検証2: 3次元の運動方程式を数値積分（糸の長さを l(t)=1-0.5 t/Tp でゆっくり縮める）
gv=1.0
def run(Tp):
    th0=np.pi/3; l0=1.0
    om0=np.sqrt(gv/(l0*np.cos(th0)))
    # 位置 (x,y,z) z 下向き正，O原点。拘束 |r|=l(t) をラグランジュ乗数で
    def lt(t): return l0-0.5*min(t,Tp)/Tp
    def ldot(t): return -0.5/Tp if t<Tp else 0.0
    def f(t,Y):
        x,y,z,vx,vy,vz=Y; rr=np.array([x,y,z]); vv=np.array([vx,vy,vz])
        lv=lt(t); ld=ldot(t)
        # 拘束 r·r=l^2 -> r·a + v·v = l ld' + ld^2 (ldd=0)
        gvec=np.array([0,0,gv])
        lam=(np.dot(vv,vv)-ld**2+np.dot(rr,gvec))/lv**2  # a = g - lam r
        a=gvec-lam*rr
        return [vx,vy,vz,*a]
    Y0=[l0*np.sin(th0),0,l0*np.cos(th0),0,om0*l0*np.sin(th0),0]
    s=solve_ivp(f,[0,Tp+30],Y0,rtol=1e-10,atol=1e-12,max_step=0.02)
    x,y,z=s.y[0],s.y[1],s.y[2]; mask=s.t>Tp+5
    rr=np.sqrt(x**2+y**2+z**2); ang=np.arctan2(np.sqrt(x**2+y**2),z)
    return np.degrees(ang[mask].min()),np.degrees(ang[mask].max())
for Tp in [300.0,1000.0]:
    print('Tp',Tp,'theta range after',run(Tp))
# 小角: theta ∝ l^{-3/4}
print('small angle exponent check:', sp.series(sp.sin(th)**4/sp.cos(th),th,0,6))
