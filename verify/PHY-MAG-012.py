# T20132-1 : 四重極磁場レンズの組（y で集束・z で発散）と2枚組の焦点
# レンズ内部はローレンツ力の3次元運動方程式を数値積分，レンズ外は直線運動（解析）
import numpy as np, sympy as sp
from scipy.integrate import solve_ivp
q,v,b,d=1.0,1.0,50.0,0.02     # m=1 で f=mv/(qbd)=1
def through_lens(s,xc,g,mass):
    # s=[x,y,z,vx,vy,vz] を x=xc-d/2 まで直進させ，レンズ内を積分して x=xc+d/2 の状態を返す
    s=np.array(s,float); t=(xc-d/2-s[0])/s[3]; s[:3]+=s[3:]*t
    def rhs(t,u):
        B=g*np.array([0,u[2],u[1]])          # (Bx,By,Bz)=(0, g z, g y)
        return np.concatenate([u[3:],q*np.cross(u[3:],B)/mass])
    ev=lambda t,u:u[0]-(xc+d/2); ev.terminal=True
    sol=solve_ivp(rhs,[0,10],s,events=ev,rtol=1e-11,atol=1e-13,max_step=d/200)
    return sol.y_events[0][0]
def cross_after(s,idx):
    # 直進して成分 idx が0になる x
    t=-s[idx]/s[3+idx]; return s[0]+s[3]*t if t>0 else None
def beam(y0,z0,lenses,mass=1.0):
    s=[-0.5,y0,z0,v,0,0]
    for xc,g in lenses: s=through_lens(s,xc,g,mass)
    return s
f=1.0/(q*b*d)
s=beam(1e-3,0,[(0,b)]); print('single y-cross x=',cross_after(s,1),' th f=',f)
s=beam(0,1e-3,[(0,b)]); print('single z virtual focus x=',s[0]-s[2]*s[3]/s[5],' th',-f)
for L in [0.5,0.3]:
    sy=beam(1e-3,0,[(0,b),(L,-b)]); sz=beam(0,1e-3,[(0,b),(L,-b)])
    print('L=',L,' y-focus behind lens2',cross_after(sy,1)-L,' th',f*(f-L)/L,' | z-focus',cross_after(sz,2)-L,' th',f*(f+L)/L)
for mm in [0.4,0.8]:
    fp=mm*f; sy=beam(1e-3,0,[(0,b),(0.5,-b)],mass=mm); c=cross_after(sy,1)
    print("m'=",mm,' y-cross behind lens2:',None if c is None else c-0.5,' th',fp*(fp-0.5)/0.5 if fp>0.5 else 'none')
F_,L_=sp.symbols('f L',positive=True)
lens=lambda fv: sp.Matrix([[1,0],[-1/fv,1]]); drift=lambda Lv: sp.Matrix([[1,Lv],[0,1]])
for nm,M in [('y',lens(-F_)*drift(L_)*lens(F_)),('z',lens(F_)*drift(L_)*lens(-F_))]:
    y,th=M*sp.Matrix([1,0]); print(nm,'back focal',sp.simplify(-y/th))
print('L=f/2:', sp.simplify((F_*(F_+L_)/L_).subs(L_,F_/2)))
