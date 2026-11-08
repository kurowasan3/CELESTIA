# P27 類題：底の穴をふさぐ球（穴の部分は大気）にはたらく水の力と浮き上がる水深
import numpy as np, sympy as sp
R_,r_,h_,rho_,g_=sp.symbols('R r h rho g',positive=True)
c_=R_-sp.sqrt(R_**2-r_**2)
F_=rho_*g_*(sp.Rational(4,3)*sp.pi*R_**3-sp.pi*c_**2*(3*R_-c_)/3)-rho_*g_*h_*sp.pi*r_**2
print('F(h)=',F_)
R,r,rho,g=0.030,0.018,1000.0,9.8
c=R-np.sqrt(R**2-r**2); V=4/3*np.pi*R**3; Vc=np.pi*c**2*(3*R-c)/3; A=np.pi*r**2
m=0.040
h0=(rho*(V-Vc)-m)/(rho*A); print('c=',c,' h0=',h0,' top of ball',2*R-c)
# 直接積分：球面のうち穴の平面（z=0, 球の中心は z=R-c）より上の部分に水圧，下の部分に大気（ゲージ0）
from scipy.integrate import quad
def Fwater(h):
    zc=R-c
    thc=np.arccos(-zc/R)   # 穴の縁（z=0）の極角
    # 外向き法線の z 成分 cos(t)，力 = -p n dS，dS = 2π R^2 sin t dt，p = ρg(h - z)
    f=lambda t: -rho*g*(h-(zc+R*np.cos(t)))*np.cos(t)*2*np.pi*R**2*np.sin(t)
    return quad(f,0,thc,epsabs=1e-13,epsrel=1e-12)[0]
for h in [0.06,h0,0.10]:
    print('h=',h,' F numeric',Fwater(h),' formula',rho*g*(V-Vc)-rho*g*h*A,' mg',m*g)
assert abs(Fwater(h0)-m*g)<1e-9
dp=rho*g*(0.20-h0); print('Δp for lift at 20 cm:',dp)
