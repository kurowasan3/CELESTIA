# PHY-OSC-054 回転するシリンダー内のピストンと等温気体
import sympy as sp, numpy as np
from scipy.integrate import solve_ivp
m,L,n,R,T,x,w=sp.symbols('m L n R T x omega',positive=True)
K=n*R*T
# つり合い m w^2 x = nRT/(L-x)
w2=K/(m*x*(L-x))
print('omega^2(x)=',sp.simplify(w2))
print('d(omega^2)/dx=0 at', sp.solve(sp.diff(w2,x),x))
# 安定条件：d/dx[ K/(L-x) - m w^2 x ] >0  at equilibrium
stab=sp.simplify((sp.diff(K/(L-x),x)-m*w2).subs({}))
print('stiffness=',sp.factor(stab))   # >0 iff x>L/2
# (3) x0=2L/3 -> x2=5L/6
w1sq=w2.subs(x,2*L/3); w2sq=w2.subs(x,5*L/6)
print('w1^2=',sp.simplify(w1sq),' w2^2=',sp.simplify(w2sq))
Wgas=K*sp.log((L-2*L/3)/(L-5*L/6))
print('W_C =',sp.simplify(Wgas))
KE=sp.Rational(1,2)*m*w2sq*(5*L/6)**2
print('motor work =',sp.simplify(Wgas+KE))
# (4) x0=L/4
wj=w2.subs(x,L/4); print('jump omega^2=',sp.simplify(wj))
xf=sp.solve(sp.Eq(w2,wj),x); print('equilibria at same omega:',xf)
Ed=sp.Rational(1,2)*m*wj*((3*L/4)**2-(L/4)**2)-K*sp.log((L-L/4)/(L-3*L/4))
print('dissipated=',sp.simplify(Ed), float(sp.simplify(Ed/K)))
# ---- numerical check (units m=L=K=1) ----
mm=1.0;LL=1.0;KK=1.0
# (3) slow ramp in rotating frame with omega(t), radial eq: m(x''-w^2 x) = -K/(L-x) (with tiny damping to suppress oscillation? no: use very slow ramp)
w1=np.sqrt(4.5); w2f=np.sqrt(36/5); Tr=4000.0
def om(t): return w1+(w2f-w1)*min(t/Tr,1.0)
def dom(t): return (w2f-w1)/Tr if t<Tr else 0.0
def f(t,y):
    x,v,Wm,Wg=y; o=om(t)
    a=o*o*x-KK/(LL-x)/mm
    N=mm*(2*o*v+dom(t)*x)       # tangential wall force
    return [v,a,N*o*x, KK/(LL-x)*v]  # motor power, work on gas = p dV(-) -> F*v
y0=[2/3,0,0,0]
s=solve_ivp(f,[0,Tr+50],y0,rtol=1e-10,atol=1e-12,max_step=1.0)
x_end=s.y[0,-1]; print('num x_end',x_end,'(5/6=0.8333)')
print('num work on gas',s.y[3,-1],' ln2=',np.log(2))
print('num motor work',s.y[2,-1],' theory',np.log(2)+2.5)
# (4) jump at constant omega with damping c
wj2=16/3; c=0.05
def g(t,y):
    x,v,E=y
    a=wj2*x-KK/(LL-x)/mm-c*v/mm
    return [v,a,c*v*v]
s2=solve_ivp(g,[0,3000],[0.25+1e-6,0,0],rtol=1e-10,atol=1e-12,max_step=0.5)
print('num final x',s2.y[0,-1],' dissipated',s2.y[2,-1],' theory',4/3-np.log(3))
