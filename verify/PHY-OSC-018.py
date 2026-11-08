# P15-1: 一定角速度で回る輪のビーズ：底からのずれ出し，最大角，モーターの仕事 W = ω ΔL
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
th,m,g,R,w,c=sp.symbols('theta m g R omega c',positive=True)
U=-m*g*R*sp.cos(th)-sp.Rational(1,2)*m*R**2*w**2*sp.sin(th)**2
print('dU:',sp.factor(sp.diff(U,th)))
# 最大角: U(th_m)=U(0)
cm=sp.symbols('c_m')
eq=sp.Eq((-m*g*R*cm-sp.Rational(1,2)*m*R**2*w**2*(1-cm**2)),-m*g*R)
print('cos th_m:',sp.solve(eq,cm))
# c=cos th0=g/(R w^2) で書き直し
sub={g:c*R*w**2}
v2=sp.simplify((2*(U.subs(th,0)-U.subs(th,sp.acos(c)))/m).subs(sub)); print('v_rel^2 at th0:',sp.factor(v2), '-> in g,R:',sp.factor(sp.simplify(v2.subs(w,sp.sqrt(g/(R*c))))))
k2=sp.simplify(sp.diff(U,th,2).subs(th,sp.acos(c)).subs(sub)/(m*R**2)); print('Omega^2:',sp.factor(k2))
sm2=1-(2*c-1)**2
W=sp.simplify((m*R**2*w**2*sm2).subs(w,sp.sqrt(g/(R*c)))); print('W =',sp.factor(W))
# 数値: 慣性系で輪が及ぼす方位角方向の力の仕事を直接積分（独立検証）
gn,Rn,mn=9.8,0.5,0.2; wn=np.sqrt(2.5*gn/Rn)  # c=0.4
cn=gn/(Rn*wn**2)
def f(t,y):
    t_,p=y; return [p, np.sin(t_)*(wn**2*np.cos(t_)-gn/Rn)]
ev=lambda t,y:y[1]; ev.terminal=True; ev.direction=-1
r=solve_ivp(f,[0,50],[1e-6,0],events=ev,rtol=1e-12,atol=1e-14,dense_output=True)
thm=r.y[0,-1]; print('th_m num',np.cos(thm),'formula',2*cn-1)
ts=np.linspace(0,r.t[-1],400001); Y=r.sol(ts)
# 方位角方向の力 F_phi = m(2 ω d(R sinθ)/dt)，速度 ω R sinθ → 仕事率
P=mn*2*wn*Rn*np.cos(Y[0])*Y[1]*wn*Rn*np.sin(Y[0])
Wnum=np.trapezoid(P,ts) if hasattr(np,'trapezoid') else np.trapz(P,ts)
print('W num',Wnum,'formula 4mgR(1-c)',4*mn*gn*Rn*(1-cn))
# 最大相対速さ
vmax=np.max(np.abs(Y[1]))*Rn; print('vmax num',vmax,'formula',np.sqrt(gn*Rn*(1-cn)**2/cn))
# 微小振動の周期
th0=np.arccos(cn); r2=solve_ivp(f,[0,20],[th0+1e-4,0],rtol=1e-12,atol=1e-15,dense_output=True)
tt=np.linspace(0,20,200001); z=r2.sol(tt)[0]-th0; zc=tt[1:][(z[:-1]>0)&(z[1:]<=0)]
print('period num',np.diff(zc).mean(),'formula',2*np.pi/(wn*np.sqrt(1-cn**2)))
