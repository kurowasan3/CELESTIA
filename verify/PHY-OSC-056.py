# PHY-OSC-056 検証：動ける台上の「円弧＋水平区間＋円弧」の往復運動
import sympy as sp, numpy as np
from scipy.integrate import quad, solve_ivp
m,M,g,R,l,phi=sp.symbols('m M g R ell phi',positive=True)
mu=m*M/(m+M)
# 水平区間での相対速さ（運動量0，エネルギー保存，微小角）
u=sp.sqrt(2*m*g*R*phi**2/2/mu)
w=sp.sqrt(m*g*R/(mu*R**2))          # E = 1/2 mu (R phidot)^2 + 1/2 m g R phi^2
print('u =',sp.simplify(u),' omega^2 =',sp.simplify(w**2))
T=2*sp.pi/w+2*l/u
print('T =',sp.simplify(T), ' check u=R w phi:',sp.simplify(u-R*w*phi))
# 逆問題
T1s,T2s,p1=sp.symbols('T1 T2 phi1',positive=True)
Tf=lambda ph: 2*sp.pi/w+2*l/(R*w*ph)
sol=sp.solve([sp.Eq(Tf(p1),T1s),sp.Eq(Tf(2*p1),T2s)],[M,l],dict=True)
print('inverse:',[ {k:sp.simplify(v) for k,v in s.items()} for s in sol])
print('expected ell:',sp.simplify(2*sp.pi*R*p1*(T1s-T2s)/(2*T2s-T1s)))
print('expected M/m:',sp.simplify(1/(4*sp.pi**2*R/(g*(2*T2s-T1s)**2)-1)))
# 数値：厳密な形状で周期を積分（非線形）して微小角公式と比較
mv,Mv,gv,Rv,lv=1.0,3.0,9.8,1.0,0.5
muv=mv*Mv/(mv+Mv)
def prof(xi):
    a=abs(xi)-lv/2
    if a<=0: return 0.0,0.0
    y=Rv-np.sqrt(Rv**2-a**2); yp=a/np.sqrt(Rv**2-a**2)*np.sign(xi); return y,yp
for ph0 in [0.02,0.05,0.1]:
    h0=Rv*(1-np.cos(ph0)); xmax=lv/2+Rv*np.sin(ph0)
    def f(xi):
        y,yp=prof(xi); return 1/np.sqrt(2*mv*gv*(h0-y)/(muv+mv*yp**2))
    # 端点の特異性は quad で処理
    Tnum=2*(quad(f,-xmax,-lv/2,limit=200)[0]+quad(f,-lv/2,lv/2)[0]+quad(f,lv/2,xmax,limit=200)[0])
    Tform=float(T.subs({m:mv,M:Mv,g:gv,R:Rv,l:lv,phi:ph0}))
    print(f'phi0={ph0}: Tnum={Tnum:.5f}  Tformula={Tform:.5f}  rel={abs(Tnum/Tform-1):.2e}')
# 台の移動範囲を運動方程式の数値積分で確認（ph0=0.05）
ph0=0.05; h0=Rv*(1-np.cos(ph0)); xmax=lv/2+Rv*np.sin(ph0)
def rhs(t,s):
    X,xi,V,vxi=s
    y,yp=prof(xi)
    a=abs(xi)-lv/2
    ypp=(Rv**2/(Rv**2-a**2)**1.5) if a>0 else 0.0
    # Lagrange: 運動量保存で X を消去した1自由度：(mu+m yp^2) xi'' + m yp ypp xi'^2 + m g yp = 0
    axi=-(mv*yp*ypp*vxi**2+mv*gv*yp)/(muv+mv*yp**2)
    return [V,vxi,-mv/(mv+Mv)*axi,axi]
s=solve_ivp(rhs,[0,30],[0,-xmax,0,0],max_step=1e-3,rtol=1e-9,atol=1e-12)
X=s.y[0]; print('platform range num =',X.max()-X.min(),' formula =',mv/(mv+Mv)*(lv+2*Rv*np.sin(ph0)))
print('max rel angle on right side =',np.arcsin((s.y[1].max()-lv/2)/Rv))
